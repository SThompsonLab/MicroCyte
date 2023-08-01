//YggData is the ThompsonLab ImageJ macro for single cell analysis of IFA images

//Default image and anchor names
CH1_default = "blue";
CH2_default = "green";
CH3_default = "red";
CH4_default = "violet";
anchor_name = "dna";
types = newArray("Batch", "Single");
//Creates the dialog box object
Dialog.create("MicroCyte settings");
Dialog.addString("Anchor name:", anchor_name);
Dialog.addString("Anchor size range:", "30-400");
Dialog.addNumber("Anchor lower threshold:", 45);
Dialog.addNumber("Anchor lower circularity (0-1):", 0.75);
Dialog.addString("Image 1 name:", CH1_default);
Dialog.addString("Image 2 name:", CH2_default);
Dialog.addString("Image 3 name:", CH3_default);
Dialog.addString("Image 4 name:", CH4_default);
Dialog.addChoice("Type:", types);  
Dialog.addCheckbox("Run peri-anchor extractions", false);
Dialog.addNumber("Enlargement factor:", 3);
Dialog.addCheckbox("Run anchor-independent extractions", false);
Dialog.show();

//Parameters are pulled to run the rest of the macro
anchorName = Dialog.getString();
anchorSize = Dialog.getString();
anchorMinThreshold = Dialog.getNumber();
circularity = Dialog.getNumber();
  
wcTarget2 = Dialog.getString();
wcTarget3 = Dialog.getString();
wcTarget4 = Dialog.getString();
wcTarget5 = Dialog.getString();
  
type = Dialog.getChoice();
runPeri = Dialog.getCheckbox();
enlargeFactor = Dialog.getNumber();
runWC = Dialog.getCheckbox();
  
wcTarget1=anchorName;
wcTarget1Threshold=anchorMinThreshold;
wcTarget1Size=anchorSize;

// If the anchor-independen box is checked, a second dialog box
//is opened to get the other variables
if (runWC){
	print("Running WC");
  	Dialog.create("Anchor-independent settings");
  	Dialog.addString(wcTarget2+" size range:", "1-infinity");
  	Dialog.addNumber(wcTarget2+" lower threshold:", 15);
  	Dialog.addNumber(wcTarget2+" lower circularity (0-1):", 0);
  
  	Dialog.addString(wcTarget3+" size range:", "1-infinity");
  	Dialog.addNumber(wcTarget3+" lower threshold:", 15);
  	Dialog.addNumber(wcTarget3+" lower circularity (0-1):", 0);
  
  	Dialog.addString(wcTarget4+" size range:", "1-infinity");
  	Dialog.addNumber(wcTarget4+" lower threshold:", 15);
  	Dialog.addNumber(wcTarget4+" lower circularity (0-1):", 0);
  
  	Dialog.addString(wcTarget5+" size range:", "1-infinity");
  	Dialog.addNumber(wcTarget5+" lower threshold:", 15);
  	Dialog.addNumber(wcTarget5+" lower circularity (0-1):", 0);
 
  	Dialog.show();
  	
  	wcTarget2Size=Dialog.getString();
  	wcTarget2Threshold=Dialog.getNumber();
  	wcTarget2Circu=Dialog.getNumber();
  	
  	wcTarget3Size=Dialog.getString();
  	wcTarget3Threshold=Dialog.getNumber();
  	wcTarget3Circu=Dialog.getNumber();
  	
  	wcTarget4Size=Dialog.getString();
  	wcTarget4Threshold=Dialog.getNumber();
  	wcTarget4Circu=Dialog.getNumber();
  	
  	wcTarget5Size=Dialog.getString();
  	wcTarget5Threshold=Dialog.getNumber();
  	wcTarget5Circu=Dialog.getNumber();
}

// First, lets make sure the measurements are set for appropriate analyses. If changes are made here, they will not be included in the default imaGen() compilation

run("Set Measurements...", "area mean standard modal min centroid center perimeter feret's integrated median area_fraction display redirect=None decimal=9");

roi_image = File.openDialog("Choose a File");
input=getDirectory("current");
parent=File.getParent(input);
input=File.getParent(parent);
input=File.getParent(input);

pList=getFileList(input);
for (i=0; i < pList.length; i++){
	if(!endsWith(pList[i], ".ijm")){
		//print(pList[i]);
		dList=getFileList(input+"/"+pList[i]);
		for (j=0; j < dList.length; j++){
			pathway = input+"/"+pList[i]+"/"+dList[j]+"/PNGS";
			iList = getFileList(pathway);
			roi_image = pathway+"/"+anchorName+".png";
			//print(roi_image);
			if(File.exists(roi_image)){
					open(roi_image);
				run("8-bit");
				setAutoThreshold("Default dark");
//run("Threshold...");
				setThreshold(anchorMinThreshold, 255);
				setOption("BlackBackground", false);
				run("Convert to Mask");
//The ROIs are rounded and split
// If you are running a high magnification (>10x) DNA image, it is recommended that you comment this out to avoid Anchor_extraction image fragementation
				run("Watershed");
//The ROIs are generated
				run("Analyze Particles...", "size="+anchorSize+" add include exclude circularity="+circularity+"-1.00");
//The image is closed
				close();			
			
// Then you generate a list of the images in that directory:
				if(!File.isDirectory(pathway+"/Anchor_extraction")){
					File.makeDirectory(pathway+"/Anchor_extraction");
				};
				inputn = pathway+"/Anchor_extraction/";
				for (k=0; k < iList.length; k++){
					if (endsWith(iList[k], ".png") & !startsWith(iList[k], "overlay")){
						open(pathway+"/"+iList[k]);
						roiManager("measure");
						saveAs("Results", inputn+replace(iList[k], ".png", ".csv"));
						run("Clear Results");
						selectWindow("Results");
						run("Close");
						close();
					};
				};
			
				if(runPeri){
// Now the ROIs re-drawn and are increased by an enlargment factor
					if(!File.isDirectory(pathway+"/PeriAnchor_extraction")){
						File.makeDirectory(pathway+"/PeriAnchor_extraction");
					};
					inputp = pathway+"/PeriAnchor_extraction/";
					open(roi_image);
					run("8-bit");
					setAutoThreshold("Default dark");
					setThreshold(anchorMinThreshold, 255);
					setOption("BlackBackground", false);
					run("Convert to Mask");
//The ROIs are rounded and split
					run("Watershed");
					counts=roiManager("count");
	
					for(l=0; l<counts; l++) {
    					roiManager("Select", l);
   					 	run("Enlarge...", "enlarge="+enlargeFactor);
    					roiManager("Update");
					};
					roiManager("deselect");
					close();
// And the images are re-analyzed with slightly larger ROIs
					for (k=0; k < iList.length; k++){
						if (endsWith(iList[k], ".png") & !startsWith(iList[k], "overlay")){
							open(pathway+"/"+iList[k]);
							roiManager("measure");
							saveAs("Results", inputp+replace(iList[k], ".png", ".csv"));
							run("Clear Results");
							selectWindow("Results");
							run("Close");
							close();
						};
					};
//Then everything is closed
					selectWindow("ROI Manager");
					run("Close");
				} else {
					selectWindow("ROI Manager");
					run("Close");
				};
			
//Now that the Anchor_extraction data is collected, Ygg will collected the Anchor_extraction indepedent data
//Since the list of images will be the same, it simply iterates through each image and collects the total ROI pixel data
				if(runWC){
					if(!File.isDirectory(pathway+"/NonAnchor_extraction")){
						File.makeDirectory(pathway+"/NonAnchor_extraction");
					};
					inputc = pathway+"/NonAnchor_extraction/";
					for (k=0; k < iList.length; k++){
						if (endsWith(iList[k], ".png") & !startsWith(iList[k], "overlay")){
							open(pathway+"/"+iList[k]);
							run("8-bit");
// Change this for deviations from default
							setAutoThreshold("Default dark");
//run("Threshold...");
							if (startsWith(iList[k], wcTarget1)){
								setThreshold(anchorMinThreshold, 255);
								setOption("BlackBackground", false);
								run("Convert to Mask");
								run("Watershed");
								run("Analyze Particles...", "size="+anchorSize+" include add exclude");
							} else if (startsWith(iList[k], wcTarget2)){
								setThreshold(wcTarget2Threshold, 255);
								setOption("BlackBackground", false);
								run("Convert to Mask");
								run("Watershed");
								run("Analyze Particles...", "size="+wcTarget2Size+" add");
							} else if (startsWith(iList[k], wcTarget3)){
								setThreshold(wcTarget3Threshold, 255);
								setOption("BlackBackground", false);
								run("Convert to Mask");
								run("Watershed");
								run("Analyze Particles...", "size="+wcTarget3Size+" add");
							} else if (startsWith(iList[k], wcTarget4)){
								setThreshold(25, 255);
								setOption("BlackBackground", false);
								run("Convert to Mask");
								run("Watershed");
								run("Analyze Particles...", "size="+wcTarget4Size+" add");
							} else if (startsWith(iList[k], wcTarget5)){
								setThreshold(wcTarget5Threshold, 255);
								setOption("BlackBackground", false);
								run("Convert to Mask");
								run("Watershed");
								run("Analyze Particles...", "size="+wcTarget5Size+" add");
							}else {
								setAutoThreshold("Default dark");
								setOption("BlackBackground", false);
								run("Convert to Mask");
								run("Watershed");
								run("Analyze Particles...", "size="+wcTargetDefaultSize+" add");
							};
							close();
							open(pathway+"/"+iList[k]);
							roiManager("measure");
							saveAs("Results", inputc+replace(iList[k], ".png", ".csv"));
							run("Clear Results");
							run("Close");
							roiManager("reset")
							if (nImages>0) {
								close();
							};
						};
					};
//Everything is closed
				selectWindow("ROI Manager");
				run("Close");
				};
			}else{
				print("Anchor image not found in "+pList[i]+" - "+dList[j]);
			};
		};
	};
}

//Saves parameters of the run
input=File.getParent(input);
f=File.open(input+"bin/MicroCyte_imj_params.txt");
print("Anchor name: "+anchorName+"\n");
print("Anchor size: "+anchorSize+"\n");
print("Anchor min threshold: "+toString(anchorMinThreshold)+"\n");
print("Anchor circularity: "+toString(circularity)+"-1\n");
if (runPeri){
	print("Peri-anchor extraction ran with an enlargement factor of "+toString(enlargeFactor)+"\n");
}
if (runWC){
	print("Anchor-independent extraction ran\n");
	print(wcTarget2+" threshold set at "+toString(wcTarget2Threshold)+", a size range of "+wcTarget2Size+", and a circularity range of "+toString(wcTarget2Circu)+"-1");
	print(wcTarget3+" threshold set at "+toString(wcTarget3Threshold)+", a size range of "+wcTarget3Size+", and a circularity range of "+toString(wcTarget3Circu)+"-1");
	print(wcTarget4+" threshold set at "+toString(wcTarget4Threshold)+", a size range of "+wcTarget4Size+", and a circularity range of "+toString(wcTarget4Circu)+"-1");
	print(wcTarget5+" threshold set at "+toString(wcTarget5Threshold)+", a size range of "+wcTarget5Size+", and a circularity range of "+toString(wcTarget5Circu)+"-1");
}