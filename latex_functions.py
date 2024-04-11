from runaway_functionsv2 import Cluster
import numpy as np
def latex_text(cluster):
    if isinstance(cluster,str):
        cv2 = Cluster(cluster, version='dr2')
        cv3 = Cluster(cluster, version='dr3')


    elif isinstance(cluster,Cluster):
        cv2 = Cluster(cluster.name, version='dr2')
        cv3 = Cluster(cluster.name, version='dr3')

    ra, dec = cv2.coordinates.ra, cv2.coordinates.dec
    ra_str, dec_str = ra.to_string(format='latex')[1:-1], dec.to_string(format='latex')[1:-1]
    dist_str = str(cv2.all['Dist'])+r"\pm"+str(cv2.all['e_Dist'])
    print(rf"{cv2.name.replace('_', ' ')} is located at $\alpha = {ra_str}, \delta = {dec_str}$ at a distance of ${dist_str}$ pc.")
    print(r"\cite{Dias2021} lists "+f"{cv2.all['N']} members in the cluster, out of which {cv3.N} meet the \hyperref[topic:dias-updated]{{membership quality criteria}}.")
    if cv3.N > 30:
        n=10
        print(rf"The 10 brightest amongst these are shown in table \ref{{tab:{cv2.name}-members}}.")
    else:
        n = cv3.N
        print(rf"These {cv3.N} members are tabulated in \ref{{tab:{cv2.name}-members}}.")
    
    print(rf"The selected {cv3.N} members are used to update the kinematic parameters of the cluster and the updated parameters are shown in table \ref{{tab:{cv2.name}-kinematics}}.")

    # runaway star name
    def runaway_name(cv3):
        runaway_names = []
        for i in range(len(cv3.runaways)):
            if not isinstance(cv3.runaways[i]['HIP'], np.ma.core.MaskedConstant):
                runaway_name = "HIP " + str(cv3.runaways[i]['HIP'])
                runaway_names.append(runaway_name)
            elif not isinstance(cv3.runaways[i]['TYC2'], np.ma.core.MaskedConstant):
                runaway_name = "TYC " + str(cv3.runaways[i]['TYC2'])
                runaway_names.append(runaway_name)
            elif not isinstance(cv3.runaways[i]['Source'], np.ma.core.MaskedConstant):
                runaway_name = "Gaia DR3 " + str(cv3.runaways[i]['Source'])
                runaway_names.append(runaway_name)
        
        if len(runaway_names) == 1:
            return f"The runaway star candidate is {runaway_names[0]}."
        elif len(runaway_names) > 1:
            return f"The runaway star candidates are {', '.join(runaway_names[:-1])}, and {runaway_names[-1]}."

    cv3.latex_table_kinematics()
    cv3.latex_table_members(n)
    print(runaway_name(cv3))
