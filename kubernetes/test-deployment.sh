#!/bin/bash
set -e

echo "Starting deployment test in Minikube..."

# Make sure minikube is running and in a good state
if ! minikube status &>/dev/null; then
  echo "Starting Minikube..."
  ./kubernetes/minikube-start.sh
else
  # Check if minikube is in a good state by trying to get the IP
  if ! minikube ip &>/dev/null; then
    echo "Minikube is in an inconsistent state. Stopping and restarting..."
    minikube stop || true
    minikube delete || true
    ./kubernetes/minikube-start.sh
  else
    echo "Minikube is already running"
  fi
fi

# Set docker env to use Minikube's Docker daemon
eval $(minikube docker-env)

# Prepare database directory for test
echo "Preparing test database directory..."
mkdir -p ./src/database/db
mkdir -p /tmp/starbase-test-data
mkdir -p /tmp/starbase-temp-data
# Set wide permissions for temp data dir
sudo chmod 777 /tmp/starbase-temp-data
# Copy existing database files for test (if any)
if [ -d "./src/database/db" ] && [ "$(ls -A ./src/database/db)" ]; then
  echo "Copying existing database files to test volume..."
  sudo cp -r ./src/database/db/* /tmp/starbase-test-data/ || true
fi

# Build the Docker image
echo "Building Docker image..."
docker build -t starbase-test:latest .

# Create a temporary persistent volume for testing
cat <<EOF > kubernetes/test-pv.yaml
apiVersion: v1
kind: PersistentVolume
metadata:
  name: starbase-test-pv
  labels:
    type: local
spec:
  storageClassName: standard
  capacity:
    storage: 1Gi
  accessModes:
    - ReadWriteOnce
  hostPath:
    path: "/tmp/starbase-test-data"
---
apiVersion: v1
kind: PersistentVolumeClaim
metadata:
  name: starbase-test-pvc
spec:
  storageClassName: standard
  accessModes:
    - ReadWriteOnce
  resources:
    requests:
      storage: 1Gi
---
apiVersion: v1
kind: PersistentVolume
metadata:
  name: starbase-temp-pv
  labels:
    type: local
spec:
  storageClassName: standard
  capacity:
    storage: 100Mi
  accessModes:
    - ReadWriteOnce
  hostPath:
    path: "/tmp/starbase-temp-data"
---
apiVersion: v1
kind: PersistentVolumeClaim
metadata:
  name: starbase-temp-pvc
spec:
  storageClassName: standard
  accessModes:
    - ReadWriteOnce
  resources:
    requests:
      storage: 100Mi
EOF

# Create a temporary test deployment file
cat <<EOF > kubernetes/test-pod.yaml
apiVersion: v1
kind: Pod
metadata:
  name: starbase-test
  labels:
    app: starbase-test
spec:
  containers:
  - name: starbase
    image: starbase-test:latest
    imagePullPolicy: Never
    ports:
    - containerPort: 8000
    env:
    - name: FLASK_ENV
      value: "test"
    - name: TEST_MODE
      value: "true"
    - name: DEBUG
      value: "true"
    - name: CELERY_SCHEDULE_DIR
      value: "/tmp/celery"
    command: ["./start-script.sh"]
    volumeMounts:
    - name: starbase-test-db
      mountPath: /app/src/database/db
  volumes:
  - name: starbase-test-db
    persistentVolumeClaim:
      claimName: starbase-test-pvc
  - name: starbase-temp-data
    persistentVolumeClaim:
      claimName: starbase-temp-pvc
EOF

# Apply the test resources
echo "Creating test persistent volume and claim..."
kubectl apply -f kubernetes/test-pv.yaml

echo "Deploying test pod..."
kubectl apply -f kubernetes/test-pod.yaml

# Wait for pod to be ready (with longer timeout for database initialization)
echo "Waiting for pod to be ready..."
kubectl wait --for=condition=Ready pod/starbase-test --timeout=180s || true

# Check pod status
echo "Pod status:"
kubectl get pod starbase-test

# Display logs to debug initialization issues
echo "Initial pod logs:"
kubectl logs starbase-test

# Wait longer for application to fully initialize
echo "Waiting for application to fully initialize..."
sleep 15

# Port forward for testing
echo "Setting up port forwarding..."
kubectl port-forward pod/starbase-test 8001:8000 &
PORT_FORWARD_PID=$!

echo "Test deployment is ready at http://localhost:8001"
echo "Running basic connectivity test..."

# Wait for the service to be available
sleep 5

# Run a basic connectivity test
if curl -s http://localhost:8001/api/cache/status | grep -q "status"; then
  echo "✅ Service is responding properly"
else
  echo "❌ Service is not responding as expected"
  echo "Checking container logs for errors:"
  kubectl logs starbase-test
fi

# Provide access to the pod for debugging
echo ""
echo "For debugging, you can access the pod with:"
echo "kubectl exec -it starbase-test -- sh"
echo ""
echo "To check application status inside the pod:"
echo "curl localhost:8000/api/cache/status"
echo ""

# Capture logs
echo "Capturing logs from test pod..."
kubectl logs starbase-test > starbase-test.log

# Clean up
echo "Press Ctrl+C to stop the test and clean up resources"

cleanup() {
  echo "Cleaning up test resources..."
  kill $PORT_FORWARD_PID 2>/dev/null || true
  kubectl delete -f kubernetes/test-pod.yaml
  kubectl delete -f kubernetes/test-pv.yaml
  rm kubernetes/test-pod.yaml
  rm kubernetes/test-pv.yaml
  echo "Test resources cleaned up"
}

trap cleanup EXIT

# Keep the script running until user presses Ctrl+C
read -r -d '' _ </dev/tty 