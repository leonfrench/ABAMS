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
package ubic.BAMSandAllen.AllenDataLoaders;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Map;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.Util;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class AllenCatalogMatrices2 {
    protected static Log log = LogFactory.getLog( AllenCatalogMatrices2.class );
    public String folder;
    public String endfix;

    public AllenCatalogMatrices2() {
        folder = SetupParameters.config.getString( "abams.dataFolder" );
        endfix = ".matrix.cache";
    }

    public String getFullFilename( String name ) {
        return folder + name + endfix;
    }

    public void loadAll() throws Exception {
        loadAllenMiner();
        loadAllenNew();
        loadAllenOrignal();
    }

    public void loadAllenNew() throws Exception {
        ObjectOutputStream out;
        NewStructureExpressionInformationLoader loader = new NewStructureExpressionInformationLoader();

        log.info( "Loading 2009 data" );

        log.info( "Loading energies" );
        out = new ObjectOutputStream( new FileOutputStream( getFullFilename( "NewEnergies" ) ) );
        out.writeObject( loader.getMatrix( NewStructureExpressionInformationLoader.energyPosition ) );
        out.close();

        log.info( "Loading levels" );
        out = new ObjectOutputStream( new FileOutputStream( getFullFilename( "NewLevels" ) ) );
        out.writeObject( loader.getMatrix( NewStructureExpressionInformationLoader.levelPosition ) );
        out.close();

        log.info( "Loading densities" );
        out = new ObjectOutputStream( new FileOutputStream( getFullFilename( "NewDensity" ) ) );
        out.writeObject( loader.getMatrix( NewStructureExpressionInformationLoader.densityPosition ) );
        out.close();
    }

    public void loadAllenOrignal() throws Exception {
        ObjectOutputStream out;
        StructureExpressionInformationLoader loader = new StructureExpressionInformationLoader();

        log.info( "Loading 2007 data" );

        log.info( "Loading energies" );
        out = new ObjectOutputStream( new FileOutputStream( getFullFilename( "OriginalEnergies" ) ) );
        out.writeObject( loader.getMatrix( StructureExpressionInformationLoader.energyPosition ) );
        out.close();

        log.info( "Loading levels" );
        out = new ObjectOutputStream( new FileOutputStream( getFullFilename( "OriginalLevels" ) ) );
        out.writeObject( loader.getMatrix( StructureExpressionInformationLoader.levelPosition ) );
        out.close();

        log.info( "Loading densities" );
        out = new ObjectOutputStream( new FileOutputStream( getFullFilename( "OriginalDensity" ) ) );
        out.writeObject( loader.getMatrix( StructureExpressionInformationLoader.densityPosition ) );
        out.close();
    }

    public void loadAllenMiner() throws Exception {
        ObjectOutputStream out;
        AllenMinerLoader loader = new AllenMinerLoader();
        log.info( "Loading allen miner data" );

        log.info( "Loading total expression level" );
        out = new ObjectOutputStream( new FileOutputStream( getFullFilename( "MinerLevel" ) ) );
        out.writeObject( loader.getMatrix( AllenMinerLoader.totalExpressionInROI ) );
        out.close();

        log.info( "Loading mean expression" );
        out = new ObjectOutputStream( new FileOutputStream( getFullFilename( "MinerMean" ) ) );
        out.writeObject( loader.getMatrix( AllenMinerLoader.meanExpressionPosition ) );
        out.close();
    }

    public DoubleMatrix<String, String> getFromDisk( String name ) throws Exception {
        ObjectInputStream in = new ObjectInputStream( new FileInputStream( getFullFilename( name ) ) );
        DoubleMatrix<String, String> result = ( DenseDoubleMatrix<String, String> ) in.readObject();
        in.close();
        return result;
    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        AllenCatalogMatrices2 loader = new AllenCatalogMatrices2();
        loader.loadAll();
        // loader.loadAllenNew();
        // DoubleMatrix<String, String> matrix;
        // String matrixName = "NewLevels";
        // matrix = loader.getFromDisk( matrixName );
        // Util.writeRTable( matrixName + ".txt", matrix );
        // matrixName = "NewDensity";
        // matrix = loader.getFromDisk( matrixName );
        // Util.writeRTable( matrixName + ".txt", matrix );
    }

}
