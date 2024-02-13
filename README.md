# Getting runaways (all functions necessary included):
```python
cluster = Cluster('Ruprecht_170')
cluster.generate_tables()
theoretical_data = theoretical_isochrone(cluster,output="table",printing=False)
fs = cluster.read_table('fs')
runaways = get_runaways(cluster,fs,theoretical_data)
display(runaways)
```

