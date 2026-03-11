#!/bin/bash
set -e

IMAGE="alphagenome-explorer"
CONTAINER="alphagenome-explorer"
APP_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "Building new image..."
docker build -t "$IMAGE" "$APP_DIR"

echo "Stopping old container..."
docker stop "$CONTAINER"
docker rm "$CONTAINER"

echo "Starting new container..."
docker run -d \
  --name "$CONTAINER" \
  --restart unless-stopped \
  -p 127.0.0.1:5060:5060 \
  "$IMAGE"

echo "Done. Container status:"
docker ps --filter "name=$CONTAINER" --format "table {{.Names}}\t{{.Status}}\t{{.Ports}}"
