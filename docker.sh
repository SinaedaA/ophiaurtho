#!/bin/bash
# ophiaurtho Docker helper script

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

print_help() {
    cat << EOF
ophiaurtho Docker Helper Script

Usage: ./docker.sh [command]

Commands:
    setup       Build containers and set up databases (first time setup)
    build       Build/rebuild containers only
    db-setup    Download and set up databases
    run         Run the pipeline with default settings
    cogs        Extract all Orthogroup fasta files from the Proteinortho graph
    report      Generate report for TF of interest (set in config.yml)
    shell       Open interactive shell in snakemake container
    clean       Remove all containers and volumes (WARNING: deletes databases!)
    status      Show status of containers and volumes
    logs        Show logs from services
    help        Show this help message

Examples:
    ./docker.sh setup              # First time setup (includes database download)
    ./docker.sh run                # Run pipeline with default config
    ./docker.sh shell              # Open shell to run custom snakemake commands
    ./docker.sh logs interproscan  # View InterProScan logs

EOF
}

check_docker() {
    if ! command -v docker &> /dev/null; then
        echo -e "${RED}Error: Docker is not installed${NC}"
        echo "Please install Docker Desktop from https://www.docker.com/products/docker-desktop"
        exit 1
    fi
    
    if ! docker info &> /dev/null; then
        echo -e "${RED}Error: Docker daemon is not running${NC}"
        echo "Please start Docker Desktop"
        exit 1
    fi
}

check_compose() {
    if ! docker compose version &> /dev/null; then
        echo -e "${RED}Error: docker-compose is not available${NC}"
        echo "Please install docker-compose or update Docker Desktop"
        exit 1
    fi
}

cmd_setup() {
    echo -e "${GREEN}Starting full setup...${NC}"
    if ! ./setup_compose.sh; then
        echo -e "${RED}Error: Platform detection failed${NC}"
        exit 1
    fi
    
    echo -e "\n${YELLOW}Step 1/2: Building containers...${NC}"
    docker compose build
    
    echo -e "\n${YELLOW}Step 2/2: Setting up databases (this may take 30min-2hrs)...${NC}"
    docker compose run --rm db-setup
    
    echo -e "\n${GREEN}Setup complete! You can now run the pipeline with: ./docker.sh run${NC}"
}

cmd_build() {
    echo -e "${GREEN}Building containers...${NC}"
    if ! ./setup_compose.sh; then
        echo -e "${RED}Error: Platform detection failed${NC}"
        exit 1
    fi
    docker compose build
    echo -e "${GREEN}Build complete!${NC}"
}

cmd_db_setup() {
    echo -e "${GREEN}Setting up databases (this may take 30min-2hrs)...${NC}"
    docker compose run --rm db-setup
    echo -e "${GREEN}Database setup complete!${NC}"
}

cmd_run() {
    echo -e "${GREEN}Running ophiaurtho pipeline...${NC}"
    docker compose run --rm snakemake bash -c "
        snakemake \
            --configfile config/config.yml \
            --cores 8 \
            --use-conda \
            --rerun-incomplete
    "
}

cmd_cogs() {
    echo -e "${GREEN}Extracting all COGs to FASTA files...${NC}"
    docker compose run --rm snakemake \
        snakemake -s workflow/rules/get_all_cogs.smk \
        --configfile config/config.yml \
        --cores 8 \
        --use-conda \
        --rerun-incomplete
}

cmd_report() {
    echo -e "${GREEN}Generating report for your favorite TF...${NC}"
    docker compose run --rm snakemake \
        snakemake -s workflow/rules/make_report.smk \
        --configfile config/config.yml \
        --cores 8 \
        --use-conda \
        --rerun-incomplete
}

cmd_shell() {
    echo -e "${GREEN}Opening interactive shell...${NC}"
    docker compose run --rm snakemake bash
}

cmd_clean() {
    echo -e "${RED}WARNING: This will delete all containers, volumes, and databases!${NC}"
    read -p "Are you sure? (yes/no): " confirm
    if [ "$confirm" = "yes" ]; then
        echo -e "${YELLOW}Cleaning up...${NC}"
        docker compose down -v
        echo -e "${GREEN}Cleanup complete${NC}"
    else
        echo "Cancelled"
    fi
}

cmd_status() {
    echo -e "${GREEN}Container Status:${NC}"
    docker compose ps
    
    echo -e "\n${GREEN}Volume Status:${NC}"
    docker volume ls | grep ophiaurtho || echo "No volumes found"
    
    echo -e "\n${GREEN}Disk Usage:${NC}"
    docker system df
}

cmd_logs() {
    service="${1:-}"
    if [ -z "$service" ]; then
        docker compose logs
    else
        docker compose logs "$service"
    fi
}

# Main script
check_docker
check_compose

case "${1:-help}" in
    setup)
        cmd_setup
        ;;
    build)
        cmd_build
        ;;
    db-setup)
        cmd_db_setup
        ;;
    run)
        cmd_run
        ;;
    shell)
        cmd_shell
        ;;
    report)
        cmd_report
        ;;
    cogs)
        cmd_cogs
        ;;
    clean)
        cmd_clean
        ;;
    status)
        cmd_status
        ;;
    logs)
        cmd_logs "${2:-}"
        ;;
    help|--help|-h)
        print_help
        ;;
    *)
        echo -e "${RED}Unknown command: $1${NC}"
        print_help
        exit 1
        ;;
esac