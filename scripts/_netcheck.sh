#!/bin/bash
echo "node: $(hostname)"
timeout 15 curl -fsS -o /dev/null -w "quay.io=%{http_code}\n" https://quay.io || echo "quay UNREACHABLE"
timeout 15 curl -fsS -o /dev/null -w "github_raw=%{http_code}\n" https://raw.githubusercontent.com || echo "githubraw UNREACHABLE"
timeout 15 curl -fsS -o /dev/null -w "s3=%{http_code}\n" https://ngi-igenomes.s3.eu-west-1.amazonaws.com || echo "s3 UNREACHABLE"
