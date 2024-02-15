# Part 0: Import necessary libraries
Import the functions necessary for the code:

```python
%run runaway_functionsv2
%matplotlib qt
```

This also imports a list of young open clusters from (based on Dias+ 2021, Gaia DR2):
```python
display(cluster_list)
```


# Part 1: Obtain the cluster class object

using the `get_cluster` function from the `runaway_functions.py`, with the `cluster_name` as the input, we obtain the parameters of the cluster.

example usage:
```python
cluster_name = 'ASCC_21'
import os
from astropy.table import Table, Column
from runaway_functions import get_cluster
cluster = get_cluster(cluster_name)
display(cluster)
```

This imports all the details for the cluster.
Various parameters of the cluster can be accessed:
- Name
- Diameter
- Number of Cluster members 

# Getting runaways (all functions necessary included):
```python
cluster = Cluster('Ruprecht_170')
cluster.generate_tables()
theoretical_data = theoretical_isochrone(cluster,output="table",printing=False)
fs = cluster.read_table('fs')
runaways = get_runaways(cluster,fs,theoretical_data)
display(runaways)
```


