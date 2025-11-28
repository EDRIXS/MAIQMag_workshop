# MAIQMag EDRIXS Course

## How to Run The Code
For the tutorial, we recommend running the code entirely online via the first option, "In cloud via Binder".

`````{tab-set}
````{tab-item} In cloud via Binder

[Open on Binder][].

````

````{tab-item} Locally via Docker

Install the [docker] application on your computer.

Download and extract the [repository].

Open the Docker Desktop app and go to the Terminal tab in the bottom right corner.
Change directory via `cd` first into the download folder and then into the primary `MAIQMag_workshop-main` folder containing the `docker-compose.yml` file and execute

```console
docker compose up
```
````

````{tab-item} Locally on linux

If you don't already have git, install it and, if needed, install pixi via

```console
curl -fsSL https://pixi.sh/install.sh | bash
```
Then the tutorials can be run via
```console
git clone https://github.com/EDRIXS/MAIQMag_workshop.git
cd MAIQMag_workshop
pixi run start
```

````
`````

Or, instead of _running_ the code, you may view the code and results by
following the links below.

## Tutorials

```{toctree}
---
maxdepth: 1
glob:
---
tutorials/intro/intro.md
tutorials/atomic/atomic_model.md
tutorials/AIM/AIM.md
```
[Open on Binder]: https://mybinder.org/v2/gh/EDRIXS/MAIQMag_workshop/HEAD?urlpath=lab/tree/tutorials/
[docker]: https://www.docker.com/products/docker-desktop/
[repository]: https://github.com/EDRIXS/MAIQMag_workshop/archive/refs/heads/main.zip
