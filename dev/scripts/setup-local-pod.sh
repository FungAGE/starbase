#!/bin/bash

# STARBASE Local Single-Pod Setup Script
# This script sets up a local environment that mimics a Kubernetes single-pod deployment

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Helper functions
print_header() {
    echo -e "${BLUE}=================================================${NC}"
    echo -e "${BLUE}  STARBASE Local Single-Pod Setup${NC}"
    echo -e "${BLUE}=================================================${NC}"
}

print_step() {
    echo -e "${GREEN}[STEP]${NC} $1"
}

print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

# Check prerequisites
check_prerequisites() {
    print_step "Checking prerequisites..."
    
    # Check for Docker
    if ! command -v docker &> /dev/null; then
        print_error "Docker is not installed. Please install Docker first."
        exit 1
    fi
    
    # Check for Docker Compose
    if ! command -v docker-compose &> /dev/null && ! docker compose version &> /dev/null; then
        print_error "Docker Compose is not installed. Please install Docker Compose first."
        exit 1
    fi
    
    # Check if Docker daemon is running
    if ! docker info &> /dev/null; then
        print_error "Docker daemon is not running. Please start Docker first."
        exit 1
    fi
    
    print_success "All prerequisites met!"
}

# Get local IP address
get_local_ip() {
    # Try different methods to get local IP
    LOCAL_IP=""
    
    # Method 1: Check common network interfaces
    for iface in eth0 en0 wlan0 wifi0; do
        IP=$(ip addr show $iface 2>/dev/null | grep 'inet ' | awk '{print $2}' | cut -d/ -f1 | head -n1)
        if [[ $IP =~ ^[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
            LOCAL_IP=$IP
            break
        fi
    done
    
    # Method 2: Use hostname -I (Linux)
    if [ -z "$LOCAL_IP" ]; then
        LOCAL_IP=$(hostname -I 2>/dev/null | awk '{print $1}')
    fi
    
    # Method 3: Use default route
    if [ -z "$LOCAL_IP" ]; then
        LOCAL_IP=$(ip route get 8.8.8.8 2>/dev/null | grep -oP 'src \K\S+')
    fi
    
    # Fallback to localhost
    if [ -z "$LOCAL_IP" ] || [[ ! $LOCAL_IP =~ ^[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
        LOCAL_IP="127.0.0.1"
        print_warning "Could not detect local IP, using localhost (127.0.0.1)"
    fi
    
    echo $LOCAL_IP
}

# Create environment file
create_env_file() {
    print_step "Creating environment configuration..."
    
    # Get local IP
    LOCAL_IP=$(get_local_ip)
    print_info "Detected local IP: $LOCAL_IP"
    
    # Create .env file from template
    if [ ! -f "dev/config/env.template" ]; then
        print_error "env.template file not found. Please ensure it exists."
        exit 1
    fi
    
    # Create .env file
    cp dev/config/env.template .env
    
    # Update .env with local IP
    sed -i "s/AUTH_DOMAIN=127.0.0.1.nip.io/AUTH_DOMAIN=$LOCAL_IP.nip.io/g" .env
    sed -i "s/HOST_IP=127.0.0.1/HOST_IP=$LOCAL_IP/g" .env
    sed -i "s/POD_IP=127.0.0.1/POD_IP=$LOCAL_IP/g" .env
    sed -i "s/KUBERNETES_SERVICE_HOST=127.0.0.1/KUBERNETES_SERVICE_HOST=$LOCAL_IP/g" .env
    
    print_success "Environment file created with local IP: $LOCAL_IP"
}

# Create required directories
create_directories() {
    print_step "Creating required directories..."
    
    # Create directories that mimic Kubernetes volume mounts
    mkdir -p src/database/cache
    mkdir -p src/database/logs
    mkdir -p src/database/db
    
    # Set proper permissions
    chmod 755 src/database/cache src/database/logs
    chmod -R 755 src/database/db
    
    print_success "Required directories created"
}

# Build Docker image
build_image() {
    print_step "Building Docker image for single-pod deployment..."
    
    # Source environment variables
    if [ -f ".env" ]; then
        export $(grep -v '^#' .env | xargs)
    fi
    
    # Build the image using local.Dockerfile
    docker build -f dev/docker/local.Dockerfile -t starbase:local-pod \
        --build-arg IPSTACK_API_KEY="${IPSTACK_API_KEY:-}" \
        --build-arg MAINTENANCE_TOKEN="${MAINTENANCE_TOKEN:-}" \
        .
    
    print_success "Docker image built successfully"
}

# Validate configuration
validate_setup() {
    print_step "Validating setup..."
    
    # Check if required files exist
    local required_files=("dev/docker/local.Dockerfile" "dev/docker/docker-compose.local.yaml" ".env" "app.py")
    for file in "${required_files[@]}"; do
        if [ ! -f "$file" ]; then
            print_error "Required file missing: $file"
            exit 1
        fi
    done
    
    # Check if directories exist
    local required_dirs=("src/database/cache" "src/database/logs" "src/database/db" "src")
    for dir in "${required_dirs[@]}"; do
        if [ ! -d "$dir" ]; then
            print_error "Required directory missing: $dir"
            exit 1
        fi
    done
    
    print_success "Setup validation completed"
}

# Start the single-pod environment
start_environment() {
    print_step "Starting STARBASE single-pod environment..."
    
    # Use docker-compose to start the environment
    if command -v docker-compose &> /dev/null; then
        docker compose -f dev/docker/docker-compose.local.yaml up -d
    else
        docker compose -f dev/docker/docker-compose.local.yaml up -d
    fi
    
    print_success "Single-pod environment started"
}

# Wait for health check
wait_for_health() {
    print_step "Waiting for application to be healthy..."
    
    local max_attempts=30
    local attempt=1
    
    while [ $attempt -le $max_attempts ]; do
        if curl -f http://localhost:8080/api/cache/status &> /dev/null; then
            print_success "Application is healthy and ready!"
            return 0
        fi
        
        print_info "Attempt $attempt/$max_attempts - waiting for application..."
        sleep 5
        ((attempt++))
    done
    
    print_warning "Health check timeout. Application may still be starting."
    print_info "Check logs with: docker logs starbase-local-pod"
}

# Display connection information
show_connection_info() {
    print_step "Connection Information"
    
    LOCAL_IP=$(get_local_ip)
    
    echo -e "${GREEN}STARBASE Single-Pod Environment is Ready!${NC}"
    echo ""
    echo -e "${BLUE}Access URLs:${NC}"
    echo -e "  Local:    http://localhost:8080"
    echo -e "  Network:  http://$LOCAL_IP:8080"
    echo -e "  Domain:   http://studio.$LOCAL_IP.nip.io:8080"
    echo ""
    echo -e "${BLUE}Pod Information:${NC}"
    echo -e "  Pod Name:      starbase-local-pod"
    echo -e "  Pod Namespace: default"
    echo -e "  Pod IP:        $LOCAL_IP"
    echo -e "  Node Name:     local-node"
    echo ""
    echo -e "${BLUE}Useful Commands:${NC}"
    echo -e "  View logs:     docker logs starbase-local-pod -f"
    echo -e "  Stop pod:      docker compose -f dev/docker/docker-compose.local.yaml down"
    echo -e "  Restart pod:   docker compose -f dev/docker/docker-compose.local.yaml restart"
    echo -e "  Shell access:  docker exec -it starbase-local-pod /bin/bash"
    echo ""
    echo -e "${YELLOW}Note:${NC} This environment simulates a single Kubernetes pod for local testing."
}

# Main execution
main() {
    print_header
    
    # Parse command line arguments
    SKIP_BUILD=false
    START_ONLY=false
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            --skip-build)
                SKIP_BUILD=true
                shift
                ;;
            --start-only)
                START_ONLY=true
                shift
                ;;
            --help)
                echo "Usage: $0 [OPTIONS]"
                echo ""
                echo "Options:"
                echo "  --skip-build    Skip Docker image building"
                echo "  --start-only    Only start the environment (skip setup)"
                echo "  --help         Show this help message"
                exit 0
                ;;
            *)
                print_error "Unknown option: $1"
                exit 1
                ;;
        esac
    done
    
    if [ "$START_ONLY" = true ]; then
        start_environment
        wait_for_health
        show_connection_info
        exit 0
    fi
    
    # Run setup steps
    check_prerequisites
    create_env_file
    create_directories
    
    if [ "$SKIP_BUILD" = false ]; then
        build_image
    fi
    
    validate_setup
    start_environment
    wait_for_health
    show_connection_info
    
    print_success "Setup completed successfully!"
}

# Run main function
main "$@"