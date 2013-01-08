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
package ubic.BAMSandAllen.FocusedAnalysis;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class ConnectivityExplore {
    private static Log log = LogFactory.getLog( ConnectivityExplore.class.getName() );

    public static void main( String args[] ) throws Exception {
        AnalyzeBAMSandAllenGenes allRegions = new AnalyzeBAMSandAllenGenes();
        allRegions.readModel( SetupParameters.getDataFolder() + "NonPropigated.rdf" );
        Direction direction = Direction.INCOMING;
        DoubleMatrix<String, String> incoming = allRegions.makeConnectionMatrix( direction );
        log.info( "Incoming connections:" + Util.zSum( incoming ) );

        direction = Direction.OUTGOING;
        DoubleMatrix<String, String> outgoing = allRegions.makeConnectionMatrix( direction );
        log.info( "outgoing connections:" + Util.zSum( outgoing ) );

        int inboth = 0;
        for ( String row : outgoing.getRowNames() ) {
            for ( String col : outgoing.getColNames() ) {
                double out = outgoing.getByKeys( row, col );
                double in = incoming.getByKeys( row, col );
                if ( out == in && in == 1d ) inboth++;
            }
        }
        log.info( "Shared connections:" + inboth );
    }

}
