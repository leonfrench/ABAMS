/*******************************************************************************
 * The ABAMS project
 * 
 * Copyright (c) 2012 University of British Columbia
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *       http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/
package ubic.BAMSandAllen.BAMSDataLoaders;

import java.io.File;
import java.io.FileInputStream;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.JenaUtil;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.StructureCatalogAnalyze;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.MatrixPairs.SimpleMatrixPair;
import ubic.BAMSandAllen.adjacency.CorrelationAdjacency;
import ubic.basecode.dataStructure.CountingMap;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.math.CorrelationStats;

import cern.colt.list.DoubleArrayList;

import com.hp.hpl.jena.rdf.model.Model;
import com.hp.hpl.jena.rdf.model.ModelFactory;
import com.hp.hpl.jena.rdf.model.Property;
import com.hp.hpl.jena.rdf.model.RDFNode;
import com.hp.hpl.jena.rdf.model.Resource;
import com.hp.hpl.jena.rdf.model.Statement;

public class BAMSFromNIFDataLoader {
    protected static Log log = LogFactory.getLog( BAMSFromNIFDataLoader.class );

    Model model;

    public BAMSFromNIFDataLoader() throws Exception {
        model = ModelFactory.createDefaultModel();
        File folder = new File( SetupParameters.config.getString( "abams.NIF.BAMS.sparql" ) );
        for ( File file : folder.listFiles() ) {
            if ( !file.getName().startsWith( "sparql (" ) ) continue;
            model.read( new FileInputStream( file ), null );
            log.info( "Model size:" + model.size() + " Added:" + file.getName() );
        }
    }

    public void getStats() {
        log.info( "Model statement size:" + model.size() );
        log.info( "Connectivity statements:"
                + model.listStatements( null, NIFVocab.source, ( RDFNode ) null ).toSet().size() );

        log.info( "   CoCoMac:" + model.listStatements( null, NIFVocab.source, "CoCoMac" ).toSet().size() );
        log.info( "   BAMS:" + model.listStatements( null, NIFVocab.source, "BAMS" ).toSet().size() );
        log.info( "   BAMS rat only:" + getBAMSConnectionResources().size() );
        log
                .info( "   ConnectomeWIKI:"
                        + model.listStatements( null, NIFVocab.source, "ConnectomeWIKI" ).toSet().size() );
        log.info( "   TemporalLobe.com:"
                + model.listStatements( null, NIFVocab.source, "TemporalLobe.com" ).toSet().size() );

        log.info( "Species:" );

        CountingMap<String> species = JenaUtil.getLiteralStringCounts( model.listStatements( null, NIFVocab.species,
                ( RDFNode ) null ) );
        log.info( "    " + species.toString() );

    }

    // returns only Rat data
    public Set<Resource> getBAMSConnectionResources() {
        List<Statement> BAMSConnectionStatements = model.listStatements( null, NIFVocab.source, "BAMS" ).toList();
        Set<Resource> result = new HashSet<Resource>();
        for ( Statement s : BAMSConnectionStatements ) {
            Resource r = s.getSubject();
            String species = JenaUtil.getStringLiteral( r, NIFVocab.species );
            if ( !species.equals( "Rat" ) ) {

            } else if ( s.getSubject().getProperty( NIFVocab.projection_strength ).getString().equals(
                    "projection_strength" ) ) {
                log.info( "Strange demo statement:" + r.toString() );
            } else {
                result.add( r );
            }
        }
        return result;
    }

    public Set<String> getConnectionRegions() {
        Set<String> regions = new HashSet<String>();
        Set<String> uniqueConnections = new HashSet<String>();
        Set<Resource> BAMSResources = getBAMSConnectionResources();
        for ( Resource r : BAMSResources ) {
            String sending = r.getProperty( NIFVocab.sending_structure ).getString();
            String recieving = JenaUtil.getStringLiteral( r, NIFVocab.receiving_structure );

            regions.add( recieving );
            regions.add( sending );
            uniqueConnections.add( recieving + sending );
        }
        log.info( "Number of regions:" + regions.size() );
        log.info( "Number of unique connections:" + uniqueConnections.size() );
        return regions;
    }

    public void printStrengthStats() {
        CountingMap<String> strenghts = new CountingMap<String>();
        for ( Resource r : getBAMSConnectionResources() ) {
            String strength = JenaUtil.getStringLiteral( r, NIFVocab.projection_strength );
            strenghts.increment( strength );
        }
        log.info( strenghts.toString() );
    }

    public DoubleMatrix<String, String> getBAMSMatrix( Direction direction, boolean useOldPartive, boolean upPropigate,
            boolean skipFibers ) throws Exception {
        if ( direction == Direction.APPENDED || useOldPartive == true ) {
            throw new RuntimeException( "Appended not supported" );
        }

        BAMSDataLoader oldBAMS = new BAMSDataLoader();
        Set<String> oldRegions = oldBAMS.getRegions();
        Set<String> regions = getConnectionRegions();

        log.info( "Intersect:" + Util.intersectSize( oldRegions, regions ) );
        for ( String inter : ( Collection<String> ) Util.subtract( regions, oldRegions ) ) {
            // log.info( "Missing:" + inter );
        }

        BAMSRDFRegionLoader XMLLoader = new BAMSRDFRegionLoader();
        Set<String> XMLRegions = XMLLoader.getRegions();

        log.info( "XML Intersect:" + Util.intersectSize( XMLRegions, regions ) );
        log.info( "XML Intersect old :" + Util.intersectSize( XMLRegions, oldRegions ) );
        log.info( "XML Subtract:" + Util.subtract( regions, XMLRegions ).size() );
        log.info( "XML Subtract:" + Util.subtract( regions, XMLRegions ) );

        Set<String> matrixRegions = new HashSet<String>();
        matrixRegions.addAll( regions );
        matrixRegions.addAll( XMLRegions );

        DoubleMatrix<String, String> connectionMatrix = new DenseDoubleMatrix<String, String>( matrixRegions.size(),
                matrixRegions.size() );
        connectionMatrix.setRowNames( new LinkedList<String>( matrixRegions ) );
        connectionMatrix.setColumnNames( new LinkedList<String>( matrixRegions ) );

        int skipped = 0;
        int alreadyMarked = 0;
        for ( Resource r : getBAMSConnectionResources() ) {
            String strength = JenaUtil.getStringLiteral( r, NIFVocab.projection_strength );
            String sending = JenaUtil.getStringLiteral( r, NIFVocab.sending_structure );
            String recieving = JenaUtil.getStringLiteral( r, NIFVocab.receiving_structure );

            if ( skipFibers && strength.equals( "fibers of passage" ) ) {// || strength.equals(
                // "light" )||
                // strength.equals( "very
                // light" ) ) {
                skipped++;
                continue;
            }

            if ( direction.equals( Direction.ANYDIRECTION ) || direction.equals( Direction.INCOMING ) ) {
                double currentValue = connectionMatrix.getByKeys( sending, recieving );
                alreadyMarked += currentValue;
                connectionMatrix.setByKeys( sending, recieving, 1d );
            }
            if ( direction.equals( Direction.ANYDIRECTION ) || direction.equals( Direction.OUTGOING ) ) {
                double currentValue = connectionMatrix.getByKeys( recieving, sending );
                alreadyMarked += currentValue;
                connectionMatrix.setByKeys( recieving, sending, 1d );
            }
        }

        log.info( "Skipped connections:" + skipped );
        log.info( "all ready marked connections:" + alreadyMarked );
        // up propigation?
        return connectionMatrix;
    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        Direction direction = Direction.OUTGOING;
        BAMSFromNIFDataLoader connectionLoader = new BAMSFromNIFDataLoader();
        connectionLoader.getStats();
        boolean skipFibers = true;
        DoubleMatrix<String, String> dataMatrix = connectionLoader.getBAMSMatrix( direction, false, false, skipFibers );
        ABAMSDataMatrix matrixA = new ABAMSDataMatrix( dataMatrix, "NIFConnectivity", new CorrelationAdjacency(
                dataMatrix ) );
        
        log.info( "new:" + matrixA.getDimensionString() );
        matrixA = matrixA.removeZeroColumns();
        matrixA = matrixA.removeZeroRows();
        log.info( "new:" + matrixA.getDimensionString() );

        skipFibers = false;
        DoubleMatrix<String, String> dataMatrixFibers = connectionLoader.getBAMSMatrix( direction, false, false,
                skipFibers );
        ABAMSDataMatrix matrixAFibers = new ABAMSDataMatrix( dataMatrixFibers, "NIFConnectivityFibers",
                new CorrelationAdjacency( dataMatrix ) );

        StructureCatalogAnalyze forMatrix = new StructureCatalogAnalyze( new BrainRegionClassSelector() );
        forMatrix.readModel( SetupParameters.getDataFolder() + "NonPropigated.rdf" );
        dataMatrix = forMatrix.makeConnectionMatrix( direction );
        ABAMSDataMatrix matrixB = new ABAMSDataMatrix( dataMatrix, "OLDConnectivity", new CorrelationAdjacency(
                dataMatrix ) );
        log.info( "old:" + matrixB.getDimensionString() );
        matrixB = matrixB.removeZeroColumns();
        matrixB = matrixB.removeZeroRows();
        log.info( "old:" + matrixB.getDimensionString() );

        SimpleMatrixPair pair2 = new SimpleMatrixPair( matrixA, matrixAFibers );
        log.info( pair2.run() );

        SimpleMatrixPair pair = new SimpleMatrixPair( matrixA, matrixB );
        log.info( pair.run() );

        int maxDiff = 0;
        String maxName = "";
        // matrixA = pair.getMatrixA();
        // matrixB = pair.getMatrixB();

        // for ( String name : matrixA.getColNames() ) {
        // double newS = Util.sum( matrixA.getColByName( name ) );
        // double oldS = Util.sum( matrixB.getColByName( name ) );
        // double diff = oldS - newS;
        // if ( diff > maxDiff ) {
        // maxName = name;
        // maxDiff = ( int ) diff;
        // log.info( "Max difference:" + maxName + " (" + maxDiff + ")" );
        // log.info( "old:" + oldS );
        // log.info( "new:" + newS );
        // }
        // }
        // log.info( "Max difference:" + maxName + " (" + maxDiff + ")" );
        log.info( "Connections OLD =" + Util.zSum( matrixB ) );
        log.info( "Connections New =" + Util.zSum( matrixA ) );

        int inNewNotOld = 0;
        for ( String colname : matrixA.getColNames() ) {
            for ( String rowname : matrixA.getRowNames() ) {
                double newValue = matrixA.getByKeys( rowname, colname );
                double oldValue = 0;
                try {
                    oldValue = matrixB.getByKeys( rowname, colname );
                } catch ( Exception e ) {
                }
                if ( newValue == 1 && oldValue == 0 ) {
                    inNewNotOld++;
                    log.info( "New not old:" + colname + "->" + rowname );
                }
            }
        }

        int intOldNotNew = 0;
        int fibersInNewConinOld = 0;
        Set<String> colnames = ( Set<String> ) Util.union( matrixB.getColNames(), matrixA.getColNames(), matrixAFibers
                .getColNames() );
        Set<String> rownames = ( Set<String> ) Util.union( matrixB.getRowNames(), matrixA.getRowNames(), matrixAFibers
                .getRowNames() );

        for ( String colname : colnames ) {
            for ( String rowname : rownames ) {
                double oldValue = 0;
                double newValue = 0;
                double newFibersSkipped = 0;
                try {
                    oldValue = matrixB.getByKeys( rowname, colname );
                    newValue = matrixA.getByKeys( rowname, colname );
                    newFibersSkipped = matrixAFibers.getByKeys( rowname, colname );
                } catch ( Exception e ) {
                }
                if ( newValue == 0 && oldValue == 1 ) {
                    intOldNotNew++;
                }
                if ( newFibersSkipped == 1 && newValue == 0 && oldValue == 1 ) {
                    fibersInNewConinOld++;
                    // log.info( colname + "->" + rowname );
                }

            }
        }
        log.info( "In old but not new:" + intOldNotNew );
        log.info( "In new but not old:" + inNewNotOld );
        log.info( "In new as just fiber and in old:" + fibersInNewConinOld );
    }
}

// http://connectivity.neuinfo.org#sending_structure
// receiving_structure
// projection_strength
// :species

class NIFVocab {
    private static Model m;
    static {
        m = ModelFactory.createDefaultModel();
    }
    protected static final String BASE = "http://connectivity.neuinfo.org#";
    public static final Property sending_structure = m.createProperty( BASE + "sending_structure" );
    public static final Property receiving_structure = m.createProperty( BASE + "receiving_structure" );
    public static final Property projection_strength = m.createProperty( BASE + "projection_strength" );
    public static final Property source = m.createProperty( BASE + "source" );
    public static final Property species = m.createProperty( BASE + "species" );
    public static final Property notes = m.createProperty( BASE + "notes" );
    public static final Property reference = m.createProperty( BASE + "reference" );

}
