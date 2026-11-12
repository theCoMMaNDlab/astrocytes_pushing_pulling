# -*- coding: utf-8 -*-
from odbAccess import openOdb
import csv
from itertools import chain

## Displacement at nodes
def Disp_Time(odb_filename,csv_filename,instance_name,step_name,node_set_name,target_times):
    odb = openOdb(path=odb_filename)
    csv_path = csv_filename 
    target_times = target_times
    if instance_name in odb.rootAssembly.instances:
        instance = odb.rootAssembly.instances[instance_name]
        
        if node_set_name in instance.nodeSets:
            node_set = instance.nodeSets[node_set_name]
            # Step 1: Extract node labels and coordinates
            node_labels = [node.label for node in node_set.nodes]
#sorted(node.label for node in node_set.nodes)
            coords_dict = {node.label: node.coordinates for node in instance.nodes}
        else:
            print("Node set {} not found in instance {}.".format(node_set_name,instance_name))
            odb.close()
            exit()
    else:
        print("Instance {} not found in ODB.".format(instance_name))
        odb.close()
        exit()
    print("Extracting node coordinates")
    # Step 2: Initialize node_data structure
    # Initialize storage
    fields_to_extract = {
        "U": ["U1", "U2", "U3"],

    }

    node_data = {
        label: {"coords": coords_dict[label], "data": []}
        for label in node_labels
    }


    for target_time in target_times:
        step = odb.steps[step_name]

        closest_frame = None
        min_time_diff = float("inf")

        for frame in step.frames:
            time_diff = abs(frame.frameValue - target_time)
            if time_diff < min_time_diff:
                min_time_diff = time_diff
                closest_frame = frame

        if closest_frame:
            print("Extracting data for time: {}".format(closest_frame.frameValue))
           # Initialize value dictionary for this time
            values_this_time = {label: None for label in node_labels}

            # Extract displacements
            disp_output = closest_frame.fieldOutputs["U"]

            for val in disp_output.values:
                if val.nodeLabel in values_this_time:
                    values_this_time[val.nodeLabel] = val.data  # U1, U2, U3

            for label in node_labels:
                node_data[label]["data"].extend(values_this_time[label])


    with open(csv_path, mode="w") as file:
        writer = csv.writer(file)
        # Build header
        base_fields = fields_to_extract["U"] #+ fields_to_extract["S"]
        header = ["Node", "X", "Y", "Z"]
        for t in target_times:
            for f in base_fields:
                header.append("{} @ {}".format(f, t))
        writer.writerow(header)



        for label in node_labels:
            flat_values = list(chain.from_iterable(
                v if isinstance(v, (list, tuple)) else [v]
                for v in node_data[label]["data"]
            ))
            writer.writerow([label] + list(node_data[label]["coords"]) + flat_values)

    print("Data saved to {}".format(csv_path))
    odb.close()
########################################################################################
if __name__ == "__main__":
    gammahat = [1,2,3]
    for ghat in gammahat:
        instance_name = "PART-1" 
        step_name = "growth"
        odb_filename = "'QuarterEllipse_pull_gammahat{}.odb".format(ghat)
        target_times = [0.0,0.1,0.15,0.2
              ,0.25,0.3,0.35,
              0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,
              0.5,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,
              0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,
              0.7,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,
              0.8,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,
              0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99]
        
        ## Cortex data
        element_set_name = "CORTEX"


        top_node_set_name = "CORTEX_TOP"

        top_disp_filename = "cortex_top_coords_disp_gammahat{}.csv".format(ghat)  
        Disp_Time(odb_filename,top_disp_filename,instance_name,step_name,top_node_set_name,target_times)
        print("Written cortex top coords & displacement file:{}".format(top_disp_filename))


