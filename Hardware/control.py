
########################
## Code snippets ONLY for the motor control and data collection
# collection of XRF and TXM data
########################


XradiaAxisNameList = ['x','y','z','t','zpx','zpy','zpz',
                      'cdx','cdy','cdz','cdpitch','cdyaw',
                      'prx','pry','prz','phx','phy','phz',
                      'energy', 'detx','dety','detz']

def collectXRF(self, filename):
    filename = filename +'.mat'
    microscope = xradia_helper.Microscope(motorwait = True, verbose = False, timer=False)
    microscope.OpenXrayShutter()

    options = getOptions(sys.argv[1:])
    print("Starting xrf system...")
    init('mercury_1.0.ini')
    start_system()
    start_run(0)
    print("waiting for 30 secs to finish the spectrum")
    time.sleep(30)
    data = get_spectrum(0)
    stop_run(0)
    print('DATA maximum: ' + str(max(data)))
    print("saving...")
    scipy.io.savemat(options.outputfile, dict(data=data))
    print("Finished xrf: ", options.outputfile)
    # plt.plot(data)
    # plt.show()
    microscope.Shutdown()

def collectTXM(self, filename):
    filename = filename +'.xrm'
    Binning = 2 # can be changed
    Exptime = 0.5 # can be changed

    microscope = xradia_helper.Microscope(motorwait = True, verbose = False, timer=False)
    if not microscope.StartGISingleAcquisition(str(os.path.join(self.filepath, filename)),
                                                          Binning,
                                                          Exptime):
                    error_message = ( 'ERROR - Could not collect image.')
                    logging.error(error_message)
                    print error_message
                    return
    print("Finished TXM: ", str(os.path.join(self.filepath, filename)))
    microscope.Shutdown()

def moveMotorAndCapture(self, filename, axisname, position):
    savefile = filename
    # self.moveBy(axisname,position)
    for index in range(len(position)):
        self.moveTo(axisname[index], position[index])
        # time.sleep(10)
    self.collectXRF(savefile)
    self.collectTXM(savefile)

Controller.moveMotorAndCapture(filename, axisname, [int(phxi), int(phyi),int(cdxi), int(cdyi)]) #not necessary to be super accurate

print 'Finished the measure, time: '+  str(endtime-starttime)

    