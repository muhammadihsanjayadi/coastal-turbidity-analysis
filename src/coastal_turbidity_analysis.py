import os
import arcpy
from arcpy.sa import Con, ExtractByMask, Float, IsNull, Raster, SetNull


# Configure the workspace and analysis inputs
project_folder = r"D:\GIS_Projects\Morowali_Turbidity_Change_Assessment\Data"
arcpy.env.workspace = project_folder
arcpy.env.overwriteOutput = True
arcpy.env.addOutputsToMap = False
arcpy.CheckOutExtension("Spatial")

study_area_path = os.path.join(project_folder, "study_area.shp")
results_folder = os.path.join(project_folder, "Turbidity_Results")

if not arcpy.Exists(study_area_path):
    raise FileNotFoundError(f"Study area does not exist: {study_area_path}")

os.makedirs(results_folder, exist_ok=True)

earlier_date = "20160402"
later_date = "20210406"
earlier_output_folder = os.path.join(results_folder, earlier_date)
later_output_folder = os.path.join(results_folder, later_date)
comparison_output_folder = os.path.join(
    results_folder, f"Comparison_{earlier_date}_to_{later_date}"
)

for output_folder in [
    earlier_output_folder,
    later_output_folder,
    comparison_output_folder,
]:
    os.makedirs(output_folder, exist_ok=True)


# Build Sentinel-2 band paths
def build_scene_paths(scene_id):
    return {
        "blue": os.path.join(project_folder, f"{scene_id}_B02_10m.jp2"),
        "green": os.path.join(project_folder, f"{scene_id}_B03_10m.jp2"),
        "red": os.path.join(project_folder, f"{scene_id}_B04_10m.jp2"),
        "nir": os.path.join(project_folder, f"{scene_id}_B08_10m.jp2"),
        "scl": os.path.join(project_folder, f"{scene_id}_SCL_20m.jp2"),
    }


earlier_scene = build_scene_paths("T51MVS_20160402T021342")
later_scene = build_scene_paths("T51MVS_20210406T021341")


# Validate Sentinel-2 inputs
def validate_scene_inputs(date_label, scene_paths):
    display_names = {
        "blue": "B02 (blue)",
        "green": "B03 (green)",
        "red": "B04 (red)",
        "nir": "B08 (NIR)",
        "scl": "SCL (20 m)",
    }
    missing_bands = [
        display_names[band_name]
        for band_name, path in scene_paths.items()
        if not arcpy.Exists(path)
    ]
    if missing_bands:
        raise FileNotFoundError(
            f"Missing Sentinel-2 input band(s) for {date_label}: "
            + ", ".join(missing_bands)
        )


# Align the Scene Classification Layer with the spectral bands
def resample_scene_classification(
    scene_classification_path, reference_band_path, output_path
):
    reference_description = arcpy.Describe(reference_band_path)
    target_cell_size = (
        f"{reference_description.meanCellWidth} "
        f"{reference_description.meanCellHeight}"
    )
    arcpy.management.Resample(
        in_raster=scene_classification_path,
        out_raster=output_path,
        cell_size=target_cell_size,
        resampling_type="NEAREST",
    )
    return Raster(output_path)


# Mask invalid and obscured pixels
def mask_invalid_pixels(spectral_bands, scene_classification):
    invalid_observations = (
        (scene_classification == 0)
        | (scene_classification == 1)
        | (scene_classification == 2)
        | (scene_classification == 3)
        | (scene_classification == 8)
        | (scene_classification == 9)
        | (scene_classification == 10)
        | (scene_classification == 11)
    )
    valid_pixel_mask = SetNull(invalid_observations, 1)
    masked_bands = {
        band_name: ExtractByMask(band, valid_pixel_mask)
        for band_name, band in spectral_bands.items()
    }
    return valid_pixel_mask, masked_bands


# Identify water with NDWI
def calculate_ndwi(green_band, nir_band, ndwi_threshold=0.00):
    green_values = Float(green_band)
    nir_values = Float(nir_band)
    ndwi_denominator = green_values + nir_values
    ndwi = SetNull(
        ndwi_denominator == 0,
        (green_values - nir_values) / ndwi_denominator,
    )
    water_mask = SetNull(ndwi <= ndwi_threshold, 1)
    return {"ndwi": ndwi, "water_mask": water_mask}


# Calculate NDTI for identified water pixels
def calculate_ndti(green_band, red_band, water_mask):
    water_green = ExtractByMask(green_band, water_mask)
    water_red = ExtractByMask(red_band, water_mask)
    green_values = Float(water_green)
    red_values = Float(water_red)
    ndti_denominator = red_values + green_values
    ndti = SetNull(
        ndti_denominator == 0,
        (red_values - green_values) / ndti_denominator,
    )
    return ndti


# Process one Sentinel-2 observation date
def process_observation_date(date_label, scene_paths, date_output_folder, study_area):
    print(f"Processing Sentinel-2 imagery for {date_label}.")
    validate_scene_inputs(date_label, scene_paths)

    arcpy.env.snapRaster = scene_paths["green"]
    arcpy.env.cellSize = scene_paths["green"]
    arcpy.env.outputCoordinateSystem = arcpy.Describe(
        scene_paths["green"]
    ).spatialReference

    original_bands = {
        "blue": Raster(scene_paths["blue"]),
        "green": Raster(scene_paths["green"]),
        "red": Raster(scene_paths["red"]),
        "nir": Raster(scene_paths["nir"]),
    }

    aligned_scl_path = os.path.join(date_output_folder, f"scl_10m_{date_label}.tif")
    aligned_scl = resample_scene_classification(
        scene_paths["scl"], scene_paths["green"], aligned_scl_path
    )

    clipped_bands = {
        band_name: ExtractByMask(band, study_area)
        for band_name, band in original_bands.items()
    }
    clipped_scl = ExtractByMask(aligned_scl, study_area)
    valid_pixel_mask, masked_bands = mask_invalid_pixels(clipped_bands, clipped_scl)

    valid_pixel_mask_path = os.path.join(
        date_output_folder, f"valid_pixels_{date_label}.tif"
    )
    valid_pixel_mask.save(valid_pixel_mask_path)

    processed_band_paths = {
        band_name: os.path.join(date_output_folder, f"{band_name}_clip_{date_label}.tif")
        for band_name in ["blue", "green", "red", "nir"]
    }
    for band_name, band in masked_bands.items():
        band.save(processed_band_paths[band_name])

    ndwi_results = calculate_ndwi(masked_bands["green"], masked_bands["nir"])
    water_ndti = calculate_ndti(
        masked_bands["green"], masked_bands["red"], ndwi_results["water_mask"]
    )

    ndwi_path = os.path.join(date_output_folder, f"ndwi_{date_label}.tif")
    water_mask_path = os.path.join(date_output_folder, f"water_mask_{date_label}.tif")
    water_ndti_path = os.path.join(date_output_folder, f"ndti_water_{date_label}.tif")

    ndwi_results["ndwi"].save(ndwi_path)
    ndwi_results["water_mask"].save(water_mask_path)
    water_ndti.save(water_ndti_path)

    true_colour_path = os.path.join(date_output_folder, f"true_color_{date_label}.tif")
    arcpy.management.CompositeBands(
        [clipped_bands["red"], clipped_bands["green"], clipped_bands["blue"]],
        true_colour_path,
    )
    print(f"True colour composite saved for {date_label}.")

    false_colour_path = os.path.join(
        date_output_folder, f"false_color_{date_label}.tif"
    )
    arcpy.management.CompositeBands(
        [clipped_bands["nir"], clipped_bands["red"], clipped_bands["green"]],
        false_colour_path,
    )
    print(f"False colour composite saved for {date_label}.")
    print(f"Finished processing {date_label}.")

    return {
        "water_mask": water_mask_path,
        "ndti_water": water_ndti_path,
        "true_color": true_colour_path,
        "false_color": false_colour_path,
    }


# Calculate change in persistent water
def calculate_persistent_water_change(
    earlier_water_mask, later_water_mask, earlier_ndti, later_ndti
):
    persistent_water_mask = Con(
        (~IsNull(earlier_water_mask)) & (~IsNull(later_water_mask)), 1
    )
    earlier_ndti_persistent = ExtractByMask(earlier_ndti, persistent_water_mask)
    later_ndti_persistent = ExtractByMask(later_ndti, persistent_water_mask)
    ndti_change = later_ndti_persistent - earlier_ndti_persistent
    return {
        "persistent_water_mask": persistent_water_mask,
        "earlier_ndti": earlier_ndti_persistent,
        "later_ndti": later_ndti_persistent,
        "ndti_change": ndti_change,
    }


# Classify NDTI change
def classify_ndti_change(ndti_change):
    return Con(
        ndti_change <= -0.05,
        1,
        Con(
            ndti_change < -0.01,
            2,
            Con(ndti_change <= 0.01, 3, Con(ndti_change < 0.05, 4, 5)),
        ),
    )


# Add descriptive labels to the NDTI change classes
def add_change_class_labels(classified_raster_path):
    class_labels = [
        "NDTI Decrease ≥ 0.05",
        "NDTI Decrease > 0.01 and < 0.05",
        "Minimal NDTI Change (±0.01)",
        "NDTI Increase > 0.01 and < 0.05",
        "NDTI Increase ≥ 0.05",
    ]
    arcpy.management.BuildRasterAttributeTable(classified_raster_path, "Overwrite")

    class_label_field = "ClassName"
    existing_fields = [field.name for field in arcpy.ListFields(classified_raster_path)]
    if class_label_field not in existing_fields:
        arcpy.management.AddField(
            classified_raster_path, class_label_field, "TEXT", field_length=50
        )

    with arcpy.da.UpdateCursor(
        classified_raster_path, ["VALUE", class_label_field]
    ) as cursor:
        for row in cursor:
            row[1] = class_labels[row[0] - 1]
            cursor.updateRow(row)


# Add outputs to ArcGIS Pro
def add_outputs_to_active_map(output_paths):
    try:
        arcgis_project = arcpy.mp.ArcGISProject("CURRENT")
        active_map = arcgis_project.activeMap
        if active_map is None:
            return

        for output in output_paths:
            if not arcpy.Exists(output):
                continue

            normalized_output_path = os.path.normcase(os.path.abspath(output))
            for layer in list(active_map.listLayers()):
                try:
                    same_data_source = (
                        layer.supports("DATASOURCE")
                        and os.path.normcase(os.path.abspath(layer.dataSource))
                        == normalized_output_path
                    )
                    if same_data_source:
                        active_map.removeLayer(layer)
                except Exception:
                    pass

            active_map.addDataFromPath(output)
    except Exception as error:
        print(f"Could not add outputs to ArcGIS Pro: {error}")


# Run the multi-date NDTI workflow
def main():
    earlier_results = process_observation_date(
        earlier_date, earlier_scene, earlier_output_folder, study_area_path
    )
    later_results = process_observation_date(
        later_date, later_scene, later_output_folder, study_area_path
    )

    # Use the earlier water-only NDTI raster as the reference grid for comparison
    arcpy.env.snapRaster = earlier_results["ndti_water"]
    arcpy.env.cellSize = earlier_results["ndti_water"]
    arcpy.env.outputCoordinateSystem = arcpy.Describe(
        earlier_results["ndti_water"]
    ).spatialReference
    arcpy.env.extent = study_area_path
    print("Comparison setup is ready.")

    change_results = calculate_persistent_water_change(
        Raster(earlier_results["water_mask"]),
        Raster(later_results["water_mask"]),
        Raster(earlier_results["ndti_water"]),
        Raster(later_results["ndti_water"]),
    )

    persistent_water_path = os.path.join(
        comparison_output_folder, f"common_water_{earlier_date}_{later_date}.tif"
    )
    earlier_ndti_path = os.path.join(
        comparison_output_folder, f"ndti_common_{earlier_date}.tif"
    )
    later_ndti_path = os.path.join(
        comparison_output_folder, f"ndti_common_{later_date}.tif"
    )
    ndti_change_path = os.path.join(
        comparison_output_folder, f"ndti_change_{earlier_date}_to_{later_date}.tif"
    )

    change_results["persistent_water_mask"].save(persistent_water_path)
    print("Common water mask saved.")
    change_results["earlier_ndti"].save(earlier_ndti_path)
    change_results["later_ndti"].save(later_ndti_path)
    print("Common-area NDTI rasters saved.")
    change_results["ndti_change"].save(ndti_change_path)
    print("NDTI change raster saved.")

    classified_change = classify_ndti_change(change_results["ndti_change"])
    classified_change_path = os.path.join(
        comparison_output_folder,
        f"ndti_change_classes_{earlier_date}_to_{later_date}.tif",
    )
    classified_change.save(classified_change_path)
    print("NDTI change-class raster saved.")
    add_change_class_labels(classified_change_path)

    outputs_to_add = [
        study_area_path,
        earlier_results["true_color"],
        earlier_results["false_color"],
        earlier_results["ndti_water"],
        later_results["true_color"],
        later_results["false_color"],
        later_results["ndti_water"],
        persistent_water_path,
        ndti_change_path,
        classified_change_path,
    ]
    add_outputs_to_active_map(outputs_to_add)
    print("Multi-date NDTI workflow completed.")


if __name__ == "__main__":
    main()
