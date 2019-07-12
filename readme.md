# read_lif

Load image from lif files to `numpy` arrays.

The code was originally written by Mathieu Leocmach in package [colloids](https://github.com/MathieuLeocmach/colloids)

## Install

- The most convenient way would be: `pip install read_lif`
- You can also Include the file `read_lif` in your project directory

## How to Use It

```python
import read_lif

reader = read_lif.Reader('lif_file.lif')
series = reader.getSeries()
chosen = series[0]  # choose first image in the lif file
image = chosen.getFrame(T=0, channel=0)  # image is a numpy array, first time point & first channel
```

For better description, see `notebooks/tutorial.ipynb`.
