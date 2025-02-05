#!/bin/bash

# Check if all 3 arguments are met
if [ $# -ne 3 ]; then
    echo "Usage: $0 <github_repo_url> <commit_hash_or_branch> <folder_path>"
    exit 1
fi

# Assign arguments to variables
REPO_URL=$1
BRANCH_OR_COMMIT=$2
FOLDER_PATH=$3

# repo name
REPO_NAME=$(basename -s .git "$REPO_URL")
mkdir -p "$REPO_NAME"
cd "$REPO_NAME" || exit 1

# Check if the directory already has a git repo
if [ ! -d ".git" ]; then
    # If there's no .git directory, initialize the repo
    git init
    git remote add origin "$REPO_URL"
else
    # If repo already exists, just fetch updates
    echo "Repository already initialized. Fetching updates..."
    git fetch origin
fi

# Enable sparse checkout if not already set
git config core.sparseCheckout true
echo "$FOLDER_PATH" > .git/info/sparse-checkout

# Checking if it's a branch (main/master) or commit
if [[ "$BRANCH_OR_COMMIT" == "main" || "$BRANCH_OR_COMMIT" == "master" ]]; then
    git checkout "$BRANCH_OR_COMMIT"
elif [[ "$BRANCH_OR_COMMIT" =~ ^[0-9a-f]{40}$ ]]; then
    git checkout "$BRANCH_OR_COMMIT"
else
    echo "Error: The provided argument '$BRANCH_OR_COMMIT' is neither a valid commit hash nor 'main'/'master' branch."
    exit 1
fi

echo "Sparse checkout completed for '$FOLDER_PATH' from '$BRANCH_OR_COMMIT'."
