"""M1: run the frozen baseline_lcdm model and save its observables."""

import classy

from cosmology import class_params, load_config
from run_class import run, save


def main():
    config = load_config()
    params = class_params(config)
    out, meta = run(params, config)
    meta["classy_version"] = classy.__version__
    save(config["name"], out, meta)
    print(f"{config['name']}: {meta['runtime_s']:.2f} s, derived = {meta['derived']}")


if __name__ == "__main__":
    main()
