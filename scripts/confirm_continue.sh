#!/bin/bash
# confirm_continue.sh

# Usage: ./confirm_continue.sh "your custom message here"
# If user confirms (y), script exits with code 0
# If user declines (n), script exits with code 1
# If called from a parent script, you can use:
#   path/to/confirm_continue.sh "Do you want to proceed?" || exit 1

msg="$1"

while true; do
  echo
  read -p "$msg (y/n) " -n 1 -r
  echo
  if [[ $REPLY =~ ^[Yy]$ ]]; then
    exit 0
  elif [[ $REPLY =~ ^[Nn]$ ]]; then
    echo "Exiting..."
    exit 1
  else
    echo "Please enter y or n."
  fi
  echo
done

