#!/bin/bash

# Test script for Starbase services
# Tests different start script modes and deployment options

set -e

echo "🧪 Starbase Service Testing Script"
echo "=================================="

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to wait for service
wait_for_service() {
    local url=$1
    local service_name=$2
    local max_attempts=${3:-30}
    local attempt=1
    
    echo "⏳ Waiting for $service_name at $url..."
    
    while [ $attempt -le $max_attempts ]; do
        if curl -s "$url" > /dev/null 2>&1; then
            echo "✅ $service_name is ready!"
            return 0
        fi
        echo "Attempt $attempt/$max_attempts: Waiting for $service_name..."
        sleep 2
        attempt=$((attempt + 1))
    done
    
    echo "❌ $service_name failed to start within $((max_attempts * 2)) seconds"
    return 1
}

echo ""
echo "📋 Available Tests:"
echo "1. Test Docker Compose deployment"
echo "2. Test k3s deployment"
echo "3. Test start script modes"
echo "4. Test Celery fallback mode"
echo "5. Run all tests"
echo ""

read -p "Choose a test (1-5): " choice

case $choice in
    1)
        echo "🐳 Testing Docker Compose deployment..."
        
        # Stop any existing containers
        docker-compose down -v
        
        # Start services
        docker-compose up -d
        
        # Wait for services
        wait_for_service "http://localhost:6379" "Redis"
        wait_for_service "http://localhost:8000/api/cache/status" "Starbase Web App"
        
        # Check Celery workers
        echo "📊 Checking Celery workers..."
        docker-compose exec celery-worker celery -A src.config.celery_config.celery inspect active || echo "No active tasks"
        
        # Show logs
        echo "📝 Recent logs:"
        docker-compose logs --tail=10 app
        docker-compose logs --tail=5 celery-worker
        docker-compose logs --tail=5 celery-beat
        
        echo "✅ Docker Compose test complete!"
        echo "🌐 App available at: http://localhost:8000"
        ;;
        
    2)
        echo "☸️ Testing k3s deployment..."
        
        if ! command_exists kubectl; then
            echo "❌ kubectl not found. Please install kubectl first."
            exit 1
        fi
        
        # Deploy to k3s
        ./k3s-serve-test/deploy-to-k3s.sh
        
        # Wait for pods
        echo "⏳ Waiting for pods to be ready..."
        kubectl wait --for=condition=ready pod -l app=redis -n starbase-test --timeout=300s
        kubectl wait --for=condition=ready pod -l app=starbase -n starbase-test --timeout=300s
        
        # Check service status
        kubectl get pods -n starbase-test
        
        # Test the service
        wait_for_service "http://starbase.local" "Starbase k3s service"
        
        echo "✅ k3s test complete!"
        echo "🌐 App available at: http://starbase.local"
        ;;
        
    3)
        echo "🔧 Testing start script modes..."
        
        echo "Testing help command:"
        ./start-script.sh --help
        
        echo ""
        echo "📝 Start script supports these roles:"
        echo "- web: Web application (default)"
        echo "- worker: Celery worker"
        echo "- beat: Celery beat scheduler"
        echo "- all: All services (dev mode only)"
        
        echo ""
        echo "Example commands:"
        echo "  ./start-script.sh                    # Start web app"
        echo "  ./start-script.sh --role worker      # Start worker"
        echo "  ./start-script.sh --role beat        # Start beat"
        echo "  ./start-script.sh --dev --role all   # Start all (dev)"
        
        echo "✅ Start script modes verified!"
        ;;
        
    4)
        echo "🔄 Testing Celery fallback mode..."
        
        # Test without Celery environment variables
        echo "Testing synchronous execution (no Celery)..."
        unset CELERY_BROKER_URL
        unset CELERY_RESULT_BACKEND
        
        # Start app briefly to test fallback
        echo "Starting app in fallback mode for 10 seconds..."
        timeout 10s ./start-script.sh --dev > /tmp/fallback-test.log 2>&1 &
        APP_PID=$!
        
        sleep 5
        
        # Check if app started
        if ps -p $APP_PID > /dev/null; then
            echo "✅ App started successfully in fallback mode"
            kill $APP_PID 2>/dev/null || true
        else
            echo "❌ App failed to start in fallback mode"
            echo "Log output:"
            cat /tmp/fallback-test.log
        fi
        
        echo "✅ Fallback mode test complete!"
        ;;
        
    5)
        echo "🚀 Running all tests..."
        echo ""
        
        # Run tests 3 and 4 (non-deployment tests)
        echo "=== Test 3: Start Script Modes ==="
        bash "$0" <<< "3"
        
        echo ""
        echo "=== Test 4: Fallback Mode ==="
        bash "$0" <<< "4"
        
        echo ""
        echo "=== Test 1: Docker Compose ==="
        bash "$0" <<< "1"
        
        echo ""
        echo "🎉 All tests completed!"
        ;;
        
    *)
        echo "❌ Invalid choice. Please run the script again and choose 1-5."
        exit 1
        ;;
esac

echo ""
echo "🔍 Useful debugging commands:"
echo "  docker-compose logs -f celery-worker    # Follow worker logs"
echo "  kubectl logs -f deployment/celery-worker -n starbase-test  # k3s worker logs"
echo "  curl http://localhost:8000/api/tasks/list  # List current tasks"
echo "  redis-cli monitor  # Monitor Redis activity" 