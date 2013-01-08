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
package ubic.BAMSandAllen;

import static ubic.BAMSandAllen.Vocabulary.direct_part_of;
import static ubic.BAMSandAllen.Vocabulary.has_direct_source;
import static ubic.BAMSandAllen.Vocabulary.has_direct_target;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.BAMSDataLoaders.BAMSDataLoader;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.ClassSelectors.ClassSelector;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.graphics.MatrixDisplay;

import com.hp.hpl.jena.ontology.OntClass;
import com.hp.hpl.jena.ontology.OntDocumentManager;
import com.hp.hpl.jena.ontology.OntModel;
import com.hp.hpl.jena.ontology.OntModelSpec;
import com.hp.hpl.jena.rdf.model.Model;
import com.hp.hpl.jena.rdf.model.ModelFactory;
import com.hp.hpl.jena.rdf.model.Property;
import com.hp.hpl.jena.rdf.model.RDFNode;
import com.hp.hpl.jena.rdf.model.Selector;
import com.hp.hpl.jena.rdf.model.Statement;
import com.hp.hpl.jena.rdf.model.StmtIterator;

public class AnalyzeBAMSandAllenGenes {
    public enum Direction {
        INCOMING, OUTGOING, ANYDIRECTION, APPENDED
    };

    private static Log log = LogFactory.getLog( AnalyzeBAMSandAllenGenes.class.getName() );

    OntModel BAMS;
    Model allenRDF;
    BAMSDataLoader BAMSData = new BAMSDataLoader();
    ClassSelector regionSelector, regionFilter;

    public AnalyzeBAMSandAllenGenes() {
        // All regions
        this( new BrainRegionClassSelector() );
    }

    public AnalyzeBAMSandAllenGenes( ClassSelector regionFilter ) {
        try {
            BAMSData = new BAMSDataLoader();
            BAMS = BAMSData.getBAMSModel();

            addAllenSubModel();

            this.regionSelector = new BrainRegionClassSelector();
            this.regionFilter = regionFilter;
        } catch ( Exception e ) {
            e.printStackTrace();
            System.exit( 1 );
        }
    }

    public void addAllenSubModel() throws FileNotFoundException {
        allenRDF = ModelFactory.createDefaultModel();
        allenRDF.read( new FileInputStream( SetupParameters.getDataFolder() + "RDFAllenToBAMSLinks.rdf" ), "" );
        BAMS.addSubModel( allenRDF );
    }

    public void writeModel( String filename ) throws IOException {
        BAMS.writeAll( new FileOutputStream( filename ), null, null );
    }

    public void readModel( String filename ) throws IOException {
        OntDocumentManager mgr = OntDocumentManager.getInstance();
        mgr.setProcessImports( false );
        OntModelSpec s = new OntModelSpec( OntModelSpec.OWL_MEM );
        s.setDocumentManager( mgr );
        BAMS = ModelFactory.createOntologyModel( s );
        BAMS.read( new FileInputStream( filename ), Vocabulary.getBAMSURI() );
    }

    public List<String> applyFilter( ClassSelector filter ) {
        Set<OntClass> filteredRegions = JenaUtil.fitlerClasses( BAMS, filter );
        return BAMSData.convertRegionstoNames( filteredRegions );
    }

    /**
     * Returns a rectangular connection matrix with incoming rows appended to outgoing (same columns)
     * 
     * @return
     */
    public DoubleMatrix<String, String> makeAppendedMatrix() {
        DoubleMatrix<String, String> in = makeConnectionMatrix( Direction.INCOMING );
        DoubleMatrix<String, String> out = makeConnectionMatrix( Direction.OUTGOING );

        // smash them together
        Set<String> colNames = new HashSet<String>( in.getColNames() );
        colNames.addAll( out.getColNames() );

        DoubleMatrix<String, String> result = new DenseDoubleMatrix<String, String>( in.rows() + out.rows(), colNames
                .size() );
        result.setColumnNames( new LinkedList<String>( colNames ) );

        for ( String inRowName : in.getRowNames() ) {
            String resultRowName = inRowName + "(incoming)";
            result.addRowName( resultRowName );
            for ( String colName : in.getColNames() ) {
                result.setByKeys( resultRowName, colName, in.getByKeys( inRowName, colName ) );
            }
        }

        // fix, code duplication
        for ( String outRowName : out.getRowNames() ) {
            String resultRowName = outRowName + "(outgoing)";
            result.addRowName( resultRowName );
            for ( String colName : out.getColNames() ) {
                result.setByKeys( resultRowName, colName, out.getByKeys( outRowName, colName ) );
            }
        }
        return result;
    }

    public DoubleMatrix<String, String> makeConnectionMatrix( Direction direction ) {
        if ( direction == Direction.APPENDED ) {
            return makeAppendedMatrix();
        }

        // only here do we use the filter
        Set<OntClass> matrixegions = JenaUtil.fitlerClasses( BAMS, regionFilter );

        List<String> enrichedRegionNames = BAMSData.convertRegionstoNames( matrixegions );
        int connections = 0;

        DoubleMatrix<String, String> connectionMatrix = new DenseDoubleMatrix<String, String>( enrichedRegionNames
                .size(), enrichedRegionNames.size() );
        connectionMatrix.setRowNames( enrichedRegionNames );
        connectionMatrix.setColumnNames( enrichedRegionNames );
        // log.info( enrichedRegionNames );

        // get all the statements for connected and filtered regions
        Selector connectedAndEnriched = new CombinedSelector( new ConnectedSelector(), new ClassSelectorWrapper(
                regionFilter ) );
        StmtIterator iter = BAMS.listStatements( connectedAndEnriched );

        while ( iter.hasNext() ) {
            Statement stmt = iter.nextStatement(); // get next statement
            OntClass subject = ( OntClass ) ( stmt.getSubject().as( OntClass.class ) ); // get the subject
            Property predicate = stmt.getPredicate(); // get the predicate
            OntClass object = ( OntClass ) ( stmt.getObject().as( OntClass.class ) ); // get the object

            /* region on the left (rows) sends output to the top (cols) - based on VanEssen */
            /*
             * this is not the normal adjacency matrix design with rows projecting to columns. The reason is the columns
             * of matrices used in mantel test have special significance
             * 
             * by setting it to incoming you will have a traditional matrix
             */
            // if ( !object.equals( subject ) ) {
            if ( predicate.equals( has_direct_source ) ) {
                // object,source is row, subject,target is the column
                connections++;
                if ( direction.equals( Direction.ANYDIRECTION ) || direction.equals( Direction.INCOMING ) ) {
                    connectionMatrix.setByKeys( BAMSData.convertClassRegionToString( object ), BAMSData
                            .convertClassRegionToString( subject ), 1d );
                }
                if ( direction.equals( Direction.ANYDIRECTION ) || direction.equals( Direction.OUTGOING ) ) {
                    connectionMatrix.setByKeys( BAMSData.convertClassRegionToString( subject ), BAMSData
                            .convertClassRegionToString( object ), 1d );
                }
            } else if ( predicate.equals( has_direct_target ) ) {
                // subject,source is row, object,target is the column
                connections++;
                if ( direction.equals( Direction.ANYDIRECTION ) || direction.equals( Direction.INCOMING ) ) {
                    connectionMatrix.setByKeys( BAMSData.convertClassRegionToString( subject ), BAMSData
                            .convertClassRegionToString( object ), 1d );
                }
                if ( direction.equals( Direction.ANYDIRECTION ) || direction.equals( Direction.OUTGOING ) ) {
                    connectionMatrix.setByKeys( BAMSData.convertClassRegionToString( object ), BAMSData
                            .convertClassRegionToString( subject ), 1d );
                }
            }
            // }
        }
        log.info( "Connections in matrix:" + connections + " (divide by two), " + connectionMatrix.rows() + " X "
                + connectionMatrix.rows() );
        return connectionMatrix;
    }

    /*
     * Go up the superclass tree until the region is one we want, or null. If current class is the right type of region,
     * it continues
     */
    public OntClass goUpPlusOne( OntClass ontClass ) {
        OntClass superClass = ontClass;
        // assume it only has one direct source
        RDFNode node = superClass.getPropertyValue( direct_part_of );
        if ( node != null ) {
            return goUp( ( OntClass ) ( node.as( OntClass.class ) ) );
        } else {
            return null;
        }
    }

    /*
     * Go up the superclass tree until the region is one we want, or null. If current class is the right type of region,
     * then its returned
     */
    public OntClass goUp( OntClass ontClass ) {// , ClassSelector stopClass ) {
        OntClass superClass = ontClass;
        while ( superClass != null && regionSelector.test( superClass ) == false ) {
            superClass = BAMSData.getParent( superClass );
        }
        return superClass;
    }

    public void propagateConnections() {
        // run a conenections statement selector
        boolean changed = false;
        do {
            StmtIterator iter = BAMS.listStatements( new ConnectedSelector() );
            List statements = iter.toList();
            changed = false;
            log.info( "Propagating.." );
            for ( Object o : statements ) {
                Statement stmt = ( Statement ) o; // get next statement
                OntClass subject = ( OntClass ) ( stmt.getSubject().as( OntClass.class ) ); // get the subject
                Property predicate = stmt.getPredicate(); // get the predicate
                OntClass object = ( OntClass ) ( stmt.getObject().as( OntClass.class ) ); // get the object

                // propigate up both sides - in a one sided manner

                // hold the object side steady
                OntClass subjectParent = goUpPlusOne( subject );
                if ( subjectParent != null && !subjectParent.equals( object ) ) {
                    // if what we are about to make does not exist
                    if ( !subjectParent.hasProperty( predicate, object ) ) {
                        changed = true;
                        subjectParent.addProperty( predicate, object );
                    }
                }

                // hold the subject side steady
                OntClass objectParent = goUpPlusOne( object );
                if ( objectParent != null && !subject.equals( objectParent ) ) {
                    // if what we are about to make does not exist
                    if ( !subject.hasProperty( predicate, objectParent ) ) {
                        changed = true;
                        subject.addProperty( predicate, objectParent );
                    }
                }
            }
        } while ( changed == true );
    }

    public void showMappings() throws Exception {
        int count = 0;
        // go through all the classes
        for ( OntClass region : JenaUtil.fitlerClasses( BAMS, new BrainRegionClassSelector() ) ) {
            OntClass superRegion = goUp( region );
            if ( superRegion == null ) {
                // System.out.println( "->Can't find valid parent" );
                count++;

            } else {
                System.out.print( region.getLabel( null ) + "(" + region.getLocalName() + ")" );
                System.out.println( "->" + BAMSData.convertClassRegionToString( superRegion ) );
            }
            // System.out.println();
        }
        System.out.println( "Number of orphans:" + count );
        /*
         * for ( String region : convertRegionstoNames( JenaUtil.fitlerClasses( BAMS, new Major17ClassSelector() ) ) ) {
         * System.out.print( region ); } System.out.println( "------------" ); for ( String region :
         * convertRegionstoNames( JenaUtil.fitlerClasses( BAMS, new Top50ClassSelector() ) ) ) { System.out.println(
         * region ); }
         */
    }

    public static void main( String args[] ) throws Exception {
        AnalyzeBAMSandAllenGenes allRegions = new AnalyzeBAMSandAllenGenes();

        Direction direction = Direction.ANYDIRECTION;
        DoubleMatrix<String, String> connectionMatrix;
        connectionMatrix = allRegions.makeConnectionMatrix( direction );

        MatrixDisplay matDisplay = new MatrixDisplay( connectionMatrix );
        matDisplay.saveImage( SetupParameters.getDataFolder() + "connectionMatrixAll.png" );

        log.info( "rows:" + connectionMatrix.rows() );
        log.info( "cols:" + connectionMatrix.columns() );

        // allRegions.showMappings();
        // allRegions.propagateConnections();
        connectionMatrix = allRegions.makeConnectionMatrix( direction );

        matDisplay = new MatrixDisplay( connectionMatrix );
        matDisplay.saveImage( SetupParameters.getDataFolder() + "connectionMatrixAllNonPropigated.png" );

        int sum = 0;
        for ( String row : connectionMatrix.getRowNames() ) {
            for ( String col : connectionMatrix.getColNames() ) {
                sum += connectionMatrix.getByKeys( row, col );
            }
        }
        log.info( sum );

        log.info( "Writing" );
        allRegions.writeModel( SetupParameters.getDataFolder() + "NonPropigated.rdf" );
        log.info( "Reading" );
        allRegions.readModel( SetupParameters.getDataFolder() + "NonPropigated.rdf" );

        connectionMatrix = allRegions.makeConnectionMatrix( direction );
        //
        // matDisplay = new MatrixDisplay( connectionMatrix );
        // matDisplay.saveImage( SetupParameters.getDataFolder() + "connectionMatrixAllPropigated.png" );
        //
        // sum = 0;
        // for ( String row : connectionMatrix.getRowNames() ) {
        // for ( String col : connectionMatrix.getColNames() ) {
        // sum += connectionMatrix.getByKeys( row, col );
        // }
        // }
        // log.info( "after read/write " + sum );
        //
        // System.exit( 1 );
        // CSVWriter myWriter = new CSVWriter( new FileWriter( SetupParameters.getDataFolder()
        // + "connectionMatrixAllDegrees.csv" ) );
        // int noconnect = 0;
        // for ( String row : connectionMatrix.getRowNames() ) {
        // int s = 0;
        // for ( String row2 : connectionMatrix.getRowNames() ) {
        // s += connectionMatrix.getByKeys( row2, row );
        // }
        // for ( String col2 : connectionMatrix.getColNames() ) {
        // s += connectionMatrix.getByKeys( row, col2 );
        // }
        // if ( s == 0 ) noconnect++;
        // myWriter.writeNext( new String[] { row.toString(), s + "" } );
        // }
        // log.info( "noconnectsio:" + noconnect );
        // myWriter.close();

    }

    /**
     * @param args
     */

}
