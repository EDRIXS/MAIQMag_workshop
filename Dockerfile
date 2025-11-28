# ── build the Pixi environment ──────────────────────────────────────────────
FROM ghcr.io/prefix-dev/pixi:0.49.0 AS build
WORKDIR /app

COPY pixi.toml  ./
RUN pixi install       # writes env to .pixi
COPY . .
EXPOSE 8888

CMD ["pixi", "run", "jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root"]



