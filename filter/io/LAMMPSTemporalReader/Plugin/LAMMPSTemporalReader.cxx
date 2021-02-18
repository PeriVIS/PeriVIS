/*=========================================================================

Program:   Visualization Toolkit

Copyright (c) Michael Bußler

=========================================================================*/

#include "LAMMPSTemporalReader.h"

#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>
#include <vtkImplicitFunction.h>
#include <vtkSmartPointer.h> 
#include <vtkCellData.h>
#include <vtkMath.h>
#include <vtkKdTreePointLocator.h>
#include <vtkSortDataArray.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <functional>
#include <string>

#ifdef WIN32
	#include <direct.h>
#else
	#include <unistd.h>
	#include <limits>
#endif

#include "qglobal.h"
#include <pqFileDialog.h>
#include "QString"

inline std::string getPath( std::string filename )
{
	std::replace( filename.begin(), filename.end(), '\\', '/'); // replace \ with /
	return filename.substr(0, filename.find_last_of("/"));
}
std::string stripExtension(std::string filename)
{
	// strip extension
	const size_t period_idx = filename.rfind('.');
	if (std::string::npos != period_idx) {
		filename.erase(period_idx);
	}
	return filename;
}
std::string stripPath( std::string filename )
{
	// strip path
	const size_t last_slash_idx = filename.find_last_of("\\/");
	if (std::string::npos != last_slash_idx) {
		filename.erase(0, last_slash_idx + 1);
	}
	return filename;
}
std::string stripPathAndExtension( std::string filename )
{
	filename = stripPath(filename);
	filename = stripExtension(filename);
	return filename;
}
bool checkFileExists( const std::string& filename)
{
    std::string fnStr(filename);
    std::replace( fnStr.begin(), fnStr.end(), '\\', '/'); // replace \ with /
    std::ifstream ftest(fnStr.c_str(), std::ifstream::in);
    if( !ftest.good()) {
        ftest.close();
        //printf("Error finding file '%s'\n", fnStr.c_str());
        return false;
    }
    ftest.close();
    return true;
}

vtkStandardNewMacro(LAMMPSTemporalReader);

LAMMPSTemporalReader::LAMMPSTemporalReader(vtkImplicitFunction *cf)
{
	//this->FileName = NULL;
	//this->BondsFileName = NULL;
	this->SetNumberOfInputPorts(0);
	this->SetNumberOfOutputPorts(1);
	this->NumberOfTimeSteps=0;
	this->bRelativePath = false;
	this->m_initialMesh = NULL;
	this->CalculateBrokenBonds = false;
	this->CalculateBondLength = true;
	this->LoadBonds = true;
}

LAMMPSTemporalReader::~LAMMPSTemporalReader()
{
    //SetFileName(NULL);
}

int LAMMPSTemporalReader::RequestInformation(vtkInformation*,
	vtkInformationVector**,
	vtkInformationVector* outputVector)
{
	vtkInformation* outInfo = outputVector->GetInformationObject(0);

	if( this->m_FileNames.empty() ) {
		vtkDebugMacro("FileName(s) not set");
		return 0;
	}

    // parse input files
    ParseLAMMPSFile();

    if( NumberOfTimeSteps == 0 ) {
        vtkErrorMacro("No timesteps found!");
        return 0;
    }

	if( this->LoadBonds && m_BondFileNames.empty() ) 
	{
        // show another file open dialog so that the user can select a bond file
        QString path = QString::fromStdString(getPath(m_FileNames[0]));
        pqFileDialog browseBondFileDialog( NULL, NULL, 
            "Select Bonds file; cancel if none!",
            path, "Bonds file (*.LAMMPS)" );
        browseBondFileDialog.setFileMode(pqFileDialog::ExistingFile);
        if(pqFileDialog::Accepted == browseBondFileDialog.exec()) 
        { 
            QStringList files = browseBondFileDialog.getSelectedFiles();
            if(!files.empty()) {
                if( m_initialMesh) {
                    m_initialMesh = NULL;
                }

                for( int i=0; i<files.size(); i++)
                    m_BondFileNames.push_back( files[i].toStdString());
            }
        }
	}

    ParseBondsFile();

	outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), 
		&this->TimestepValues[0], this->NumberOfTimeSteps);
	
	double timeRange[2] = { 
		TimestepValues[0], 
		TimestepValues[NumberOfTimeSteps-1]
	};

	outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),
		timeRange, 2);

	return 1;
}

int LAMMPSTemporalReader::RequestData(vtkInformation *vtkNotUsed(request),
	vtkInformationVector **vtkNotUsed(inputVector),
	vtkInformationVector *outputVector)
{

	//this->DebugOn();

	// get the info object
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	vtkPolyData *output = vtkPolyData::
		SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	if (output == 0) {
		vtkErrorMacro("Error while creating reader output.");
	}

	int timestepToLoad = 0;
	if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())
		&& NumberOfTimeSteps > 0 ) {

			double requestedTimeValue = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

			timestepToLoad = FindClosestTimeStep(requestedTimeValue);
			output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), TimestepValues[timestepToLoad]);
	}

	//---------------------------------------------------------------------------
	if( !m_initialMesh )
	{
		this->m_initialMesh = vtkSmartPointer<vtkPolyData>::New();
		ReadAtoms(0, this->m_initialMesh);
		if(this->LoadBonds)
            ReadBonds(0, this->m_initialMesh);
		m_initialMesh->BuildLinks();
	}
	if( timestepToLoad == 0 )
	{
		output->ShallowCopy(this->m_initialMesh);
	}
	else
	{
		ReadAtoms(timestepToLoad, output);
        if(this->LoadBonds)
    		ReadBonds(timestepToLoad, output);
	}

    output->Squeeze();
	//---------------------------------------------------------------------------

	return 1;
}

void LAMMPSTemporalReader::ParseLAMMPSFile()
{
	using namespace std;

	m_timesteps.clear();
	TimestepValues.clear();

    for( int i=0; i<m_FileNames.size(); i++)
    {

        ifstream inFile;
        inFile.open(m_FileNames[i].c_str(), ifstream::in);

        // get length of file
        inFile.seekg(0, inFile.end);
        size_t filesize = inFile.tellg();
        inFile.seekg(0,inFile.beg);

        Timestep ts;
        ts.atomStartLine = 0;
        ts.time = -1;
        ts.nBonds = -1;
        ts.bondsStartLine = 0;
        
        string line;
        uint64_t lineId = 0;

        while (inFile.good()) 
        {
            getline(inFile, line);
            lineId++;

            vector<string> tokens;
            tokenize(line, tokens);

            if( tokens.size() == 0 ) {
                continue;
            }

            if (tokens[0] == "ITEM:" )
            {
                if( tokens.size() > 1 && tokens[1] == "TIMESTEP") 
                {
                    // save current timestep
                    //if( ts.time >= 0 )
                    //    m_timesteps.push_back(ts);

                    // read time
                    getline(inFile, line);lineId++;

                    // init new timestep
                    ts = Timestep();
                    ts.time = stoi(line);
                    ts.nBonds = -1;
                    ts.bondsStartLine = 0;
                    ts.atomStartLine = 0;
                    ts.dataFieldNames.clear();
                    ts.nAtoms = -1;
                    ts.nDataFields = 0;

                    //update progress
                    size_t pos = inFile.tellg();
                    this->UpdateProgress(0.5*(pos/(float)filesize));

                }
                else if ( tokens.size() > 3 && tokens[1] == "NUMBER" && tokens[2] == "OF" && tokens[3] == "ATOMS") 
                {
                    // read number of atoms
                    getline(inFile, line);lineId++;
                    ts.nAtoms = stoi(line);
                }
                else if ( tokens.size() > 3 && tokens[1] == "BOX" && tokens[2] == "BOUNDS")
                {
                    // read bounding box
                    string::size_type sz;
                    for( int i=0; i<3; i++){
                        getline(inFile, line);lineId++;
                        ts.box[2*i+0] = stof(line, &sz);
                        ts.box[2*i+1] = stof(line.substr(sz));
                    }
                }
                else if ( tokens.size() > 1 && tokens[1] == "ATOMS" )
                {
                    // read atom data, 
                    // first 4 entries are reserved for: id, x, y, z
                    ts.dataFieldNames.push_back("id");
                    int i=2;
                    while( tokens[++i] != "x" && i<(tokens.size()-2));
                    if( i==tokens.size()){
                        vtkErrorMacro("could not locate position of x field in dump data!");
                        return;
                    }
                    ts.posField[0] = (i-2); i++;
                    ts.posField[1] = (i-2); i++;
                    ts.posField[2] = (i-2); i++;

                    for( ; i<tokens.size(); i++)
                    {
                        ts.dataFieldNames.push_back(tokens[i]);
                    }
                    ts.nDataFields = ts.dataFieldNames.size();

                    // save line offset
                    ts.atomStartLine = lineId;

                    // skip atom data
                    //for (int currLineNumber = 0; currLineNumber < ts.nAtoms; ++currLineNumber){
                    //    inFile.ignore(numeric_limits<streamsize>::max(), inFile.widen('\n'));
                    //}
                    //lineId+=ts.nAtoms;
                    break;
                }
            }
        }

        inFile.close();

        // save last timestep
        if( ts.time >= 0 )
            m_timesteps.push_back(ts);
    }

	
	NumberOfTimeSteps = m_timesteps.size();

	for( int i=0; i<NumberOfTimeSteps; i++)
		TimestepValues.push_back(m_timesteps[i].time);

}

void LAMMPSTemporalReader::ParseBondsFile()
{
	if( m_BondFileNames.empty() ){
		this->UpdateProgress(1.0); // done!
		return;
	}

	using namespace std;

    for( int i=0; i<m_BondFileNames.size(); i++)
    {

        ifstream inFile;
        inFile.open(m_BondFileNames[i], ifstream::in);

        // get length of file
        inFile.seekg(0, inFile.end);
        size_t filesize = inFile.tellg();
        inFile.seekg(0,inFile.beg);

        string line;
        uint64_t lineId = 0;

        Timestep* ts = NULL;
        
        while (inFile.good()) 
        {
            getline(inFile, line);lineId++;

            vector<string> tokens;
            tokenize(line, tokens);

            if( tokens.size() == 0 ) {
                continue;
            }

            if (tokens[0] == "ITEM:" )
            {
                if( tokens.size() > 1 && tokens[1] == "TIMESTEP")
                {
                    // read time
                    getline(inFile, line);lineId++;
                    int timeFound = stoi(line);
                    int timeStepId = FindClosestTimeStep(timeFound);
                    if( timeFound == m_timesteps[timeStepId].time ){
                        ts = &m_timesteps[timeStepId];
                    } else {
                        ts = NULL;
                    }

                    //update progress
                    size_t pos = inFile.tellg();
                    this->UpdateProgress(0.5+(0.5*(pos/(float)filesize)));

                }
                else if ( tokens.size() > 3 && tokens[1] == "NUMBER" && tokens[2] == "OF" && tokens[3] == "ENTRIES") 
                {
                    // read number of atoms
                    getline(inFile, line);lineId++;
                    int numBonds = stoi(line);
                    if( ts ) {
                        ts->nBonds = numBonds;
                    } else {
                        // skip entries
                        for (int currLineNumber = 0; currLineNumber < numBonds; ++currLineNumber){
                            inFile.ignore(numeric_limits<streamsize>::max(), inFile.widen('\n'));
                        }
                        lineId += numBonds;
                    }
                }
                else if ( tokens.size() > 1 && tokens[1] == "ENTRIES" )
                {
                    if( ts ) {
                        ts->bondField[0]=0;
                        ts->bondField[1]=1;
                        if( tokens.size() >2 && tokens[2] == "index" ) {
                            ts->bondField[0]=1;
                            ts->bondField[1]=2;
                        }
                        ts->bondsStartLine = lineId;


                        // skip entries
                        //for (int currLineNumber = 0; currLineNumber < ts->nBonds; ++currLineNumber){
                        //    inFile.ignore(numeric_limits<streamsize>::max(), inFile.widen('\n'));
                        //}
                        //lineId += ts->nBonds;

                        break; // break file parsing after header read completely
                    }
                }
            }
        }
        inFile.close();
    }
}


void LAMMPSTemporalReader::ReadAtoms( int id, vtkPolyData* pd )
{
	using namespace std;

	ifstream inFile;
    std::string FileName = m_FileNames[id];
	inFile.open(FileName, ifstream::in);

	Timestep ts = m_timesteps[id];
	if( ts.atomStartLine == 0)
		return;

	// atom positions
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(ts.nAtoms);
	
	
	vtkSmartPointer<vtkCellArray> newVerts = NULL;
	if(!this->LoadBonds) 
		newVerts = vtkSmartPointer<vtkCellArray>::New();

	// init atom data arrays
	std::vector<vtkSmartPointer<vtkFloatArray> > arrs;
	for( int dataId = 0; dataId< ts.nDataFields; dataId++)
	{
		vtkSmartPointer<vtkFloatArray> data = vtkSmartPointer<vtkFloatArray>::New();
		data->SetName(ts.dataFieldNames[dataId].c_str());
		data->SetNumberOfComponents(1);
		data->SetNumberOfTuples(ts.nAtoms);
		arrs.push_back(data);
	}

	// skip to start line of data
	for (uint64_t currLineNumber = 0; currLineNumber < ts.atomStartLine; ++currLineNumber){
		inFile.ignore(numeric_limits<streamsize>::max(), inFile.widen('\n'));
	}

	string line;
	vtkIdType lineId = 0;
	float values[20];
	while (inFile.good() && lineId < ts.nAtoms ) {
			
		getline(inFile, line);
		stringstream iss(line);

		int id=-1;
		float tmp,x,y,z;
		iss >> id;
		values[0] = id;

		for( int i=1; i<ts.posField[0]; i++ )
			iss >> tmp;

		iss >> x >> y >> z;
		for( int i=1; i<ts.nDataFields; i++)
			iss >> values[i];
		iss >> ws; // skip whitespace at the end

		// check for reading errors
		if( iss.fail() ) {
			vtkDebugMacro("Error while reading line \""<<line<<"\" ("<<lineId<<"): not enough values!");
		} else if( !iss.eof() ) {
			vtkDebugMacro("Error while reading line \""<<line<<"\" ("<<lineId<<"): too many entries!");
		} else if( id<1 || id > ts.nAtoms ) {
			vtkDebugMacro("Error while reading line \""<<line<<"\" ("<<lineId<<"): Wrong particle id!");
			continue;
		} else {
			// everything was read properly
		}

		vtkIdType ptId = id-1;
		points->SetPoint(ptId, x, y, z );

		if( !this->LoadBonds )
			newVerts->InsertNextCell(1, &ptId );

		for( int i=0; i<ts.nDataFields; i++)
			arrs[i]->SetValue( ptId, values[i]);
		
		this->UpdateProgress ((lineId)/(float)(ts.nAtoms+(this->LoadBonds)?ts.nBonds:0));

		lineId++;
	}

	inFile.close();

	pd->SetPoints(points);
	if( !this->LoadBonds ){
		pd->SetVerts(newVerts);
	}
	foreach( vtkSmartPointer<vtkFloatArray> a, arrs ) {
		pd->GetPointData()->AddArray(a);
	}
}

void LAMMPSTemporalReader::ReadBonds(int timestepToLoad, vtkPolyData * pd)
{
	const Timestep& ts = m_timesteps[timestepToLoad];

	if( timestepToLoad >= m_BondFileNames.size() || ts.nBonds < 0 || ts.bondsStartLine == 0)
		return;

	if( this->LoadBonds ) {
		pd->SetVerts(NULL);
	}

	using namespace std;

	ifstream inFile;
	inFile.open(m_BondFileNames[timestepToLoad], ifstream::in);

	// bonds
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

    // bond active flag: 1:active; 0:broken
	vtkSmartPointer<vtkFloatArray> bondActive = vtkSmartPointer<vtkFloatArray>::New();;
	bondActive->SetName("Bond Active");
	bondActive->SetNumberOfComponents(1);
    
    // length of bond in current timestep
    vtkSmartPointer<vtkFloatArray> bondLength = vtkSmartPointer<vtkFloatArray>::New();
    bondLength->SetName("Bond Length");
    bondLength->SetNumberOfComponents(1);
    
    // length of bond in reference timestep
    vtkSmartPointer<vtkFloatArray> refBondLength = vtkSmartPointer<vtkFloatArray>::New();
    refBondLength->SetName("Reference Bond Length");
    refBondLength->SetNumberOfComponents(1);


	// skip to start line of data
	for (uint64_t currLineNumber = 0; currLineNumber < ts.bondsStartLine; ++currLineNumber){
		inFile.ignore(numeric_limits<streamsize>::max(), inFile.widen('\n'));
	}

	string line;
	int numLines = ts.nBonds;

	while (inFile.good() && numLines>0 ) 
	{
		getline(inFile, line);
		stringstream iss(line);
        int data[3] = {0,0,0};
        for( size_t i=0; i<=ts.bondField[1]; i++)
		    iss >> data[i];

        int count  =data[0];
        int a1     =data[ts.bondField[0]]; 
        int a2     =data[ts.bondField[1]];

		if( iss.fail() ){
			std::cout<< "Failed to read bonds from: " << line << " (line: " << ts.nBonds-numLines+1 << ")" << std::endl;
			numLines--;
			continue;
		}

		if( a1>0 && a2>0 )
		{
			vtkIdType bondIndex = lines->InsertNextCell(2);
			lines->InsertCellPoint(a1-1);
			lines->InsertCellPoint(a2-1);
		} else {
			std::cout<< "Illegal Bond found: " << line << " (line: " << ts.nBonds-numLines+1 << ")"<< std::endl;
		}
			
		numLines--;
		this->UpdateProgress ((ts.nAtoms+(ts.nBonds-numLines))/(float)(ts.nAtoms+ts.nBonds));

	} // end while

	// done with input file
	inFile.close();

	pd->SetLines(lines);

	// alternativ method to calculate broken bonds, requires consistent particle ids
	if( this->CalculateBrokenBonds ) 
	{
		int numBondsRef = m_initialMesh->GetNumberOfLines();
		bondActive->SetNumberOfTuples(numBondsRef);
		
		if( timestepToLoad == 0 ) 
		{
			for( int cellId=0; cellId<numBondsRef; cellId++ )
				bondActive->SetValue(cellId, 1.0); // at T=0 all bonds are active (1)
		} 
		else 
		{

			// initalize result array
			for( int cellId=0; cellId<numBondsRef; cellId++ ) {
				bondActive->SetValue(cellId, 0.0); // initially all bonds are non-active (1)
            }

			// use reference lines set with changes to current state
			vtkSmartPointer<vtkCellArray> reflines = vtkSmartPointer<vtkCellArray>::New();
			reflines->DeepCopy(this->m_initialMesh->GetLines()); // use bonds of initial timestep
			pd->BuildLinks(); // current state

			// bonds per point id in current and ref data set
			vtkSmartPointer<vtkIdList> cellIdsRef = vtkSmartPointer<vtkIdList>::New();
			vtkSmartPointer<vtkIdList> cellIdsCurr = vtkSmartPointer<vtkIdList>::New();

			// now we need to compare the vertices of each bond in current and ref config against each other
			// assume that there are no more particles and bonds in the current config than in ref

			int numPoints = pd->GetNumberOfPoints();
			for( vtkIdType ptId=0; ptId<numPoints; ptId++) {

				vtkSmartPointer<vtkIdList> bondsIds = vtkSmartPointer<vtkIdList>::New();
				vtkSmartPointer<vtkIdList> neighborsCurr = getConnectedPoints(pd, ptId);
				vtkSmartPointer<vtkIdList> neighborsRef = getConnectedPoints(m_initialMesh, ptId, bondsIds);

				// compare ids to ref config
				for(vtkIdType i = 0; i<neighborsRef->GetNumberOfIds(); i++)
				{
					vtkIdType neighIdRef = neighborsRef->GetId(i);

					// check if id is also in current neighbor set
					for( vtkIdType j = 0; j<neighborsCurr->GetNumberOfIds(); j++)
					{
						vtkIdType neighIdCurr = neighborsCurr->GetId(j);
						if( neighIdCurr == neighIdRef ) { // bond is still active
							vtkIdType bondId = bondsIds->GetId(i);
							bondActive->SetValue(bondId, 1.0);
							break;
						}
					}
				}
				this->UpdateProgress (0.75+0.25*ptId/(float)numPoints);
			}
			pd->SetLines(reflines);
        }
        pd->GetCellData()->AddArray(bondActive);
	}

	if(this->CalculateBondLength) 
	{
		pd->BuildCells();
		size_t numBonds = pd->GetNumberOfCells();

		bondLength->SetNumberOfTuples(numBonds);

		for(vtkIdType bondId = 0; bondId < numBonds; bondId++)
		{
			vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
			pd->GetCellPoints( bondId, pointIdList);
			double p0[3], p1[3];
            /* now use point coordinates from current state */
			pd->GetPoint( pointIdList->GetId(0), p0);
			pd->GetPoint( pointIdList->GetId(1), p1);
			double d2 = vtkMath::Distance2BetweenPoints(p0,p1);
			bondLength->SetValue(bondId, sqrt(d2));
		}

		pd->GetCellData()->AddArray(bondLength);
	}


    if( this->CalculateRefBondLength )
    {
        pd->BuildCells();
        size_t numBonds = pd->GetNumberOfCells();

        refBondLength->SetNumberOfTuples(numBonds);

        for( int bondId=0; bondId<numBonds; bondId++ ) 
        {
            vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
            pd->GetCellPoints( bondId, pointIdList);
            double p0[3], p1[3];
            /* now use point coordinates from reference state */
            m_initialMesh->GetPoint( pointIdList->GetId(0), p0);
            m_initialMesh->GetPoint( pointIdList->GetId(1), p1);
            double d2 = vtkMath::Distance2BetweenPoints(p0,p1);
            refBondLength->SetValue(bondId, sqrt(d2));
        }

        pd->GetCellData()->AddArray(refBondLength);
    }

    if( timestepToLoad > 0 && this->FilterBondsByLength && this->CalculateRefBondLength )
    {
        // new bonds, filtered by length
        vtkSmartPointer<vtkCellArray> newLines = vtkSmartPointer<vtkCellArray>::New();
        
        // new cell data arrays
        vtkSmartPointer<vtkFloatArray> bondActive_new = vtkSmartPointer<vtkFloatArray>::New();;
        bondActive_new->SetName("Bond Active");
        bondActive_new->SetNumberOfComponents(1);
        vtkSmartPointer<vtkFloatArray> bondLength_new = vtkSmartPointer<vtkFloatArray>::New();
        bondLength_new->SetName("Bond Length");
        bondLength_new->SetNumberOfComponents(1);
        vtkSmartPointer<vtkFloatArray> refBondLength_new = vtkSmartPointer<vtkFloatArray>::New();
        refBondLength_new->SetName("Reference Bond Length");
        refBondLength_new->SetNumberOfComponents(1);

        // current lines
        vtkCellArray* lines = pd->GetLines();
        size_t numBonds = lines->GetNumberOfCells();
        vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();

        // actual filtering
        for(vtkIdType bondId = 0; bondId < numBonds; bondId++)
        {
            float length = refBondLength->GetValue(bondId);
            if( length <= this->MaxBondLength )
            {
                pd->GetCellPoints(bondId, pts);
                //lines->GetCell(bondId, pts);
                newLines->InsertNextCell(pts);

                // copy cell data
                if( CalculateBrokenBonds)
                    bondActive_new->InsertNextValue( bondActive->GetValue(bondId));
                if( CalculateBondLength )
                    bondLength_new->InsertNextValue( bondLength->GetValue(bondId));
                refBondLength_new->InsertNextValue( refBondLength->GetValue(bondId));
            }
        }

        // replace data with new data
        pd->SetLines(newLines);
        if( CalculateBrokenBonds ){
            pd->GetCellData()->RemoveArray("Bond Active");
            pd->GetCellData()->AddArray(bondActive_new);
        }
        if( CalculateBondLength ){
            pd->GetCellData()->RemoveArray("Bond Length");
            pd->GetCellData()->AddArray(bondLength_new);
        }
        pd->GetCellData()->RemoveArray("Reference Bond Length");
        pd->GetCellData()->AddArray(refBondLength_new);
    }
}

vtkSmartPointer<vtkIdList> LAMMPSTemporalReader::getConnectedPoints( vtkPolyData * pd, vtkIdType ptId, vtkSmartPointer<vtkIdList> bondIds /*= NULL*/)
{
	vtkSmartPointer<vtkIdList> cellIds  = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> neighPts = vtkSmartPointer<vtkIdList>::New();

	pd->GetPointCells( ptId, cellIds);

	for(vtkIdType i = 0; i < cellIds->GetNumberOfIds(); i++)
	{
		vtkIdType cellId = cellIds->GetId(i);

		vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
		pd->GetCellPoints( cellId, pointIdList);

		vtkIdType ptIdN = pointIdList->GetId(0);
		if( ptIdN == ptId)
			ptIdN = pointIdList->GetId(1);

		neighPts->InsertNextId(ptIdN);

		if( bondIds )
			bondIds->InsertNextId(cellId);
	}

	return neighPts;
}

void LAMMPSTemporalReader::WriteHeaderFile()
{
	//// estimate header filename
 //   std::string FileName = m_FileNames[0];
	//std::string headerFilename = stripPath(FileName);
	//std::string path = getPath(FileName);

	//headerFilename = path + "/."+headerFilename;

	//ofstream outFile;
	//outFile.open(headerFilename, std::ofstream::out);

	//if( outFile.good() ){

	//	// write bonds filename
	//	std::string bondsFilename = "";
 //       if( this->BondsFileName) 
 //           stripPath(this->BondsFileName);
	//	outFile << "# " << bondsFilename << endl;

	//	// write timesteps to header file
	//	int i=0;
	//	foreach( Timestep ts, m_timesteps){

	//		std::stringstream outLine;
	//		outLine << ts.time 			 << " "
	//			<< ts.box[0]			 << " "
	//			<< ts.box[1]			 << " "
	//			<< ts.box[2]			 << " "
	//			<< ts.box[3]			 << " "
	//			<< ts.box[4]			 << " "
	//			<< ts.box[5]			 << " "
 //               << ts.nBonds			 << " "
 //               << ts.bondsStartLine	 << " "
 //               << ts.bondField[0]       << " "
 //               << ts.bondField[1]       << " "
	//			<< ts.nAtoms 			 << " "
	//			<< ts.atomStartLine		 << " "
 //               << ts.posField[0]		 << " "
 //               << ts.posField[1]		 << " "
 //               << ts.posField[2]		 << " "
 //               << ts.nDataFields        ;

	//		foreach( std::string s, ts.dataFieldNames )
	//			outLine	<< " " << s;

	//		if( ++i < m_timesteps.size() )
	//			outLine << endl;
	//		outFile << outLine.str();
	//	}
	//}
	//outFile.close();
}

bool LAMMPSTemporalReader::ReadHeaderFile()
{
	//using namespace std;

	//// estimate header filename
 //   std::string FileName = m_FileNames[0];
	//std::string headerFilename = stripPath(FileName);
	//std::string path = getPath(FileName);

	//headerFilename = path + "/."+headerFilename;

	//ifstream inFile(headerFilename);

	//if( inFile.is_open() ) {

	//	string line;

	//	while( inFile.good() ){

	//		getline(inFile, line);

 //           if( line[0] == '#') {

 //               // read bonds file name, prevent file open dialog
 //               vector<string> ts;
 //               tokenize(line,ts);
 //               SetBondsFileName(NULL);
 //               if( ts.size() > 1 ) {
 //                   // read bonds filename
 //                   string bondsFileName = path+'/'+ts[1];
 //                   if( checkFileExists(bondsFileName))
 //                       SetBondsFileName(bondsFileName.c_str());
 //               }
 //           } else {

 //               stringstream iss(line);

 //               Timestep ts;
 //               iss >> ts.time
 //                   >> ts.box[0] >> ts.box[1] >> ts.box[2] 
 //                   >> ts.box[3] >> ts.box[4] >> ts.box[5]
 //                   >> ts.nBonds
 //                   >> ts.bondsStartLine
 //                   >> ts.bondField[0] >> ts.bondField[1]
 //                   >> ts.nAtoms
 //                   >> ts.atomStartLine
 //                   >> ts.posField[0] >> ts.posField[1] >> ts.posField[2]
 //                   >> ts.nDataFields;

 //               ts.dataFieldNames.resize(ts.nDataFields);
 //               for( int i=0; i<ts.nDataFields; i++ )
 //                   iss	>> ts.dataFieldNames[i];

 //               m_timesteps.push_back(ts);
 //           }
	//	}

	//	NumberOfTimeSteps = m_timesteps.size();
	//	for( int i=0; i<NumberOfTimeSteps; i++)
	//		TimestepValues.push_back(m_timesteps[i].time);

	//	return true;
	//}

	return false;
}

int LAMMPSTemporalReader::FindClosestTimeStep(double requestedTimeValue)
{
	int ts = 0;
	double mindist = std::abs(TimestepValues[0] - requestedTimeValue);

	for (int i = 0; i < NumberOfTimeSteps; i++) {

		double tsv = TimestepValues[i];
		double dist = std::abs(tsv - requestedTimeValue);
		if (dist < mindist) {
			mindist = dist;
			ts = i;
		}
	}
	return ts;
}

void LAMMPSTemporalReader::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);

	//os << indent << "File Name: "
	//	<< (this->FileName ? this->FileName : "(none)") << "\n";
}

void LAMMPSTemporalReader::AddFileName(const char* fname)
{
    m_FileNames.push_back(fname);
}

void LAMMPSTemporalReader::RemoveAllFileNames()
{
    m_FileNames.clear();
}

void LAMMPSTemporalReader::tokenize(const std::string line, std::vector<std::string> &tokens)
{
	std::stringstream iss(line);
	std::string token;
	while( iss >> token)
		tokens.push_back(token);
}
