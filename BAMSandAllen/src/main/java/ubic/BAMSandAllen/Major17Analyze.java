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

import static ubic.BAMSandAllen.Vocabulary.hasAllen17Mapping;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.ClassSelectors.Major17ClassSelector;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.graphics.MatrixDisplay;

import com.hp.hpl.jena.ontology.OntClass;
import com.hp.hpl.jena.rdf.model.RDFNode;

public class Major17Analyze extends AnalyzeBAMSandAllenGenes {
    private static Log log = LogFactory.getLog( Major17Analyze.class.getName() );

    public Major17Analyze() {
        super( new Major17ClassSelector() );
        // regionSelector = ;
    }

    public static void main( String[] args ) throws Exception {
        Direction direction = Direction.INCOMING;

        Major17Analyze major17 = new Major17Analyze();
        major17.propagateConnections();

        DoubleMatrix<String, String> connectionMatrix;
        connectionMatrix = major17.makeConnectionMatrix(direction);
        // connectionMatrix.makePNG( "data\\connectionMatrixMajor17.png", 12 );
        MatrixDisplay matDisplay = new MatrixDisplay( connectionMatrix );
        matDisplay.saveImage( "data\\connectionMatrixMajor17.png" );

        major17.showMappings();
    }

    public String classToString( OntClass ontClass ) {
        // if its a top17 region are different
        if ( ontClass.hasProperty( hasAllen17Mapping ) ) {
            RDFNode allenLiteral = ontClass.getPropertyValue( hasAllen17Mapping );
            return allenLiteral.toString();
        } else {
            return BAMSData.convertClassRegionToString( ontClass );
        }
    }

}
