from paraview.simple import *
import os
import numpy as np
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Define the base path to the folder containing the .nek5000 file
#base_path = 'C:\\Users\\jga91\\OneDrive - Danmarks Tekniske Universitet\\SEM Speciale\\Kode\\SEM_speciale\\Jesper\\2D Unsteady Navier-stokes\\Figures to stability comp\\test'

main_directory = 'C:\\Path\\To\\Main\\Directory'

# Loop through each subdirectory in the main directory
for subdir, dirs, files in os.walk(main_directory):
    # Automatically detect the .nek5000 file in the directory
    nek5000_files = [f for f in os.listdir(base_path) if f.endswith('.nek5000')]
    nek5000_file = next((file for file in files if file.endswith('.nek5000')), None)
    if nek5000_file:
        base_path = subdir

        # Define paths for the velocity and pressure picture folders
        vel_pic_folder = os.path.join(base_path, 'vel_pic2')
        pres_pic_folder = os.path.join(base_path, 'pres_pic2')

        # Ensure the directories exist
        os.makedirs(vel_pic_folder, exist_ok=True)
        os.makedirs(pres_pic_folder, exist_ok=True)

        # Define paths for the images
        vel_pic_path = os.path.join(vel_pic_folder, 'vel.png')
        pres_pic_path = os.path.join(pres_pic_folder, 'pres.png')




        # create a new 'Nek5000 Reader'
        domain = Nek5000Reader(registrationName=nek5000_files[0], FileName=nek5000_file)

        # get the material library
        materialLibrary1 = GetMaterialLibrary()

        # get animation scene
        animationScene1 = GetAnimationScene()

        # update animation scene based on data timesteps
        animationScene1.UpdateAnimationUsingDataTimeSteps()

        # get active view
        renderView1 = GetActiveViewOrCreate('RenderView')
        # show data in view
        domain_display = Show(domain, renderView1, 'UnstructuredGridRepresentation')

        # trace defaults for the display properties.
        domain_display.Representation = 'Surface'

        # reset view to fit data
        renderView1.ResetCamera(False, 0.9)

        #changing interaction mode based on data extents
        renderView1.InteractionMode = '2D'
        renderView1.CameraPosition = [1.024999976158142, 0.20499999821186066, 6.867499840259552]
        renderView1.CameraFocalPoint = [1.024999976158142, 0.20499999821186066, 0.0]

        # get the material library
        materialLibrary1 = GetMaterialLibrary()

        # update the view to ensure updated data information
        renderView1.Update()

        # create a new 'Calculator'
        calculator1 = Calculator(registrationName='Calculator1', Input=domain)

        # Properties modified on calculator1
        calculator1.Function = 'sqrt((coordsX-0.2)^2 + (coordsY-0.2)^2)'

        # show data in view
        calculator1Display = Show(calculator1, renderView1, 'UnstructuredGridRepresentation')

        # trace defaults for the display properties.
        calculator1Display.Representation = 'Surface'

        # hide data in view
        Hide(domain, renderView1)

        # show color bar/color legend
        calculator1Display.SetScalarBarVisibility(renderView1, True)

        # update the view to ensure updated data information
        renderView1.Update()

        # get color transfer function/color map for 'Result'
        resultLUT = GetColorTransferFunction('Result')

        # get opacity transfer function/opacity map for 'Result'
        resultPWF = GetOpacityTransferFunction('Result')

        # get 2D transfer function for 'Result'
        resultTF2D = GetTransferFunction2D('Result')

        # create a new 'Threshold'
        threshold1 = Threshold(registrationName='Threshold1', Input=calculator1)

        # Properties modified on threshold1
        threshold1.UpperThreshold = 0.05

        # show data in view
        threshold1Display = Show(threshold1, renderView1, 'UnstructuredGridRepresentation')

        # trace defaults for the display properties.
        threshold1Display.Representation = 'Surface'

        # hide data in view
        Hide(calculator1, renderView1)

        # show color bar/color legend
        threshold1Display.SetScalarBarVisibility(renderView1, True)

        # update the view to ensure updated data information
        renderView1.Update()

        # turn off scalar coloring
        ColorBy(threshold1Display, None)

        # Hide the scalar bar for this color map if no visible data is colored by it.
        HideScalarBarIfNotNeeded(resultLUT, renderView1)

        # change solid color
        threshold1Display.AmbientColor = [0.0, 0.0, 0.0]
        threshold1Display.DiffuseColor = [0.0, 0.0, 0.0]

        # set active source
        SetActiveSource(domain)

        # show data in view
        domain_display = Show(domain, renderView1, 'UnstructuredGridRepresentation')

        # set scalar coloring
        ColorBy(domain_display, ('POINTS', 'Velocity', 'Magnitude'))

        # rescale color and/or opacity maps used to include current data range
        domain_display.RescaleTransferFunctionToDataRange(True, False)

        # show color bar/color legend
        domain_display.SetScalarBarVisibility(renderView1, True)

        # get color transfer function/color map for 'Velocity'
        velocityLUT = GetColorTransferFunction('Velocity')

        # get opacity transfer function/opacity map for 'Velocity'
        velocityPWF = GetOpacityTransferFunction('Velocity')

        # get 2D transfer function for 'Velocity'
        velocityTF2D = GetTransferFunction2D('Velocity')

        # get color legend/bar for velocityLUT in view renderView1
        velocityLUTColorBar = GetScalarBar(velocityLUT, renderView1)

        # change scalar bar placement
        velocityLUTColorBar.Orientation = 'Horizontal'
        velocityLUTColorBar.WindowLocation = 'Any Location'
        velocityLUTColorBar.Position = [0.3571967654986524, 0]
        velocityLUTColorBar.ScalarBarLength = 0.32999999999999907

        # set active source
        SetActiveSource(domain)


        # current camera placement for renderView1
        renderView1.InteractionMode = '2D'
        renderView1.CameraPosition = [1.061649465222839, 0.09331238189085266, 4.038725110599645]
        renderView1.CameraFocalPoint = [1.061649465222839, 0.09331238189085266, 0.0]
        renderView1.CameraParallelScale = 0.40300800585909785
        # Adjust camera

        # current camera placement for renderView1
        renderView1.InteractionMode = '2D'
        renderView1.CameraPosition = [1.061649465222839, 0.09331238189085266, 4.038725110599645]
        renderView1.CameraFocalPoint = [1.061649465222839, 0.09331238189085266, 0.0]
        renderView1.CameraParallelScale = 0.40300800585909785

        # Properties modified on velocityMagnitudeLUTColorBar
        velocityLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
        velocityLUTColorBar.TitleFontFamily = 'Times'
        velocityLUTColorBar.TitleFontSize = 25
        velocityLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
        velocityLUTColorBar.LabelFontFamily = 'Times'
        velocityLUTColorBar.LabelFontSize = 25
        # Adjust camera

        # Rescale transfer function
        velocityLUT.RescaleTransferFunction(0.0, 2.3)

        # Rescale transfer function
        velocityPWF.RescaleTransferFunction(0.0, 2.3)

        # Rescale 2D transfer function
        velocityTF2D.RescaleTransferFunction(0.0, 2.3, 0.0, 1.0)
            
        # get layout
        layout1 = GetLayout()

        # split cell
        layout1.SplitVertical(0, 0.5)
        # layout/tab size in pixels
        layout1.SetSize(1854, 267)

        layout1.SetSize(2000, 1000)

        #-----------------------------------
        # saving camera placements for views

        # current camera placement for renderView1
        renderView1.InteractionMode = '2D'
        renderView1.CameraPosition = [1.0, 0.09107615418189408, 6.6921304299024635]
        renderView1.CameraFocalPoint = [1.0, 0.09107615418189408, 0.0]
        renderView1.CameraParallelScale = 0.40300800585909774

        # save animation
        SaveAnimation(filename=vel_pic_path, viewOrLayout=renderView1, location=16, ImageResolution=[1854, 534],
            TransparentBackground=1,
            FrameWindow=[0, 150])
        # Adjust camera


        # create a new 'Contour'
        contour1 = Contour(registrationName='Contour1', Input=domain)

        # set active source
        SetActiveSource(domain)
        # Adjust camera

        # set scalar coloring
        ColorBy(domain_display, ('POINTS', 'Pressure'))

        # Hide the scalar bar for this color map if no visible data is colored by it.
        HideScalarBarIfNotNeeded(velocityLUT, renderView1)

        # rescale color and/or opacity maps used to include current data range
        domain_display.RescaleTransferFunctionToDataRange(True, False)

        # show color bar/color legend
        domain_display.SetScalarBarVisibility(renderView1, True)

        # get 2D transfer function for 'Pressure'
        pressureTF2D = GetTransferFunction2D('Pressure')

        # get color transfer function/color map for 'Pressure'
        pressureLUT = GetColorTransferFunction('Pressure')
        pressureLUT.TransferFunction2D = pressureTF2D
        pressureLUT.RGBPoints = [-0.8684472441673279, 0.231373, 0.298039, 0.752941, 0.3456191122531891, 0.865003, 0.865003, 0.865003, 1.559685468673706, 0.705882, 0.0156863, 0.14902]
        pressureLUT.ScalarRangeInitialized = 1.0

        # get opacity transfer function/opacity map for 'Pressure'
        pressurePWF = GetOpacityTransferFunction('Pressure')
        pressurePWF.Points = [-0.8684472441673279, 0.0, 0.5, 0.0, 1.559685468673706, 1.0, 0.5, 0.0]
        pressurePWF.ScalarRangeInitialized = 1
        # Adjust camera

        # set active source
        SetActiveSource(contour1)

        contour1.Isosurfaces = np.linspace(-1.4,1.7,25)
        contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

        # turn off scalar coloring
        ColorBy(contour1Display, ('POINTS', None))

        # Hide the scalar bar for this color map if no visible data is colored by it.
        HideScalarBarIfNotNeeded(pressureLUT, renderView1)
        # Adjust camera

        # Properties modified on contour1Display
        contour1Display.LineWidth = 2.0

        # get color transfer function/color map for 'Pressure'
        pressureLUT = GetColorTransferFunction('Pressure')
        pressureLUT.TransferFunction2D = pressureTF2D
        pressureLUT.RGBPoints = [-0.8684472441673279, 0.231373, 0.298039, 0.752941, 0.3456191122531891, 0.865003, 0.865003, 0.865003, 1.559685468673706, 0.705882, 0.0156863, 0.14902]
        pressureLUT.ScalarRangeInitialized = 1.0

        # get color legend/bar for pressureLUT in view renderView1
        pressureLUTColorBar = GetScalarBar(pressureLUT, renderView1)
        pressureLUTColorBar.Title = 'Pressure'
        pressureLUTColorBar.ComponentTitle = ''

        # change scalar bar placement
        pressureLUTColorBar.Orientation = 'Horizontal'
        pressureLUTColorBar.WindowLocation = 'Any Location'
        pressureLUTColorBar.Position = [0.28785329018338707, 0]
        pressureLUTColorBar.ScalarBarLength = 0.3300000000000001

        # Properties modified on velocityLUTColorBar
        pressureLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
        pressureLUTColorBar.TitleFontFamily = 'Times'
        pressureLUTColorBar.TitleFontSize = 25
        pressureLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
        pressureLUTColorBar.LabelFontFamily = 'Times'
        pressureLUTColorBar.LabelFontSize = 25



        # save animation
        SaveAnimation(filename=pres_pic_path, viewOrLayout=renderView1, location=16, ImageResolution=[1854, 534],
            TransparentBackground=1,
            FrameWindow=[0, 150])
        # Adjust camera