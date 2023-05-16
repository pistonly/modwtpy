# modwtpy
modwt in python

find the detail from the matlab Documentation:
http://www.mathworks.com/help/wavelet/ref/modwt.html,
http://www.mathworks.com/help/wavelet/ref/modwtmra.html,

and the wavelet is imported from pywt:
http://www.pybytes.com/pywavelets/ref/wavelets.html

## requirements

```shell 
pip install PyWavelets
```

## Example:
```python
from modwt import modwt, modwtmra
import pandas as pd


gdpdata = pd.read_csv('GDPcomponents.csv')
wt = modwt(gdpdata['govtexp'], 'db2', 5)
wtmra = modwtmra(wt, 'db2')
```
