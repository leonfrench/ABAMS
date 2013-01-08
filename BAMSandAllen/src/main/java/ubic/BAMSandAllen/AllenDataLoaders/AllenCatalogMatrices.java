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
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

/*
 * Basically a cache for the Allen expressio ninformation data, so I don't have to reparse the original data everytime
 */
public class AllenCatalogMatrices {
    protected static Log log = LogFactory.getLog( AllenCatalogMatrices.class );
    
    

    private DoubleMatrix<String, String> levels;
    private DoubleMatrix<String, String> energies;
    private DoubleMatrix<String, String> minerEnergies;
    private DoubleMatrix<String, String> minerMeanEnergies;
    private DoubleMatrix<String, String> densities;

    // number of unique ID's minus the number of no gene expression info
    public static final String densityFilename = SetupParameters.config.getString( "abams.dataFolder" )
            + "AllenCataLogDensities.matrix.cache";
    public static final String levelFilename = SetupParameters.config.getString( "abams.dataFolder" )
            + "AllenCatalogLevels.matrix.cache";
    public static final String energyFilename = SetupParameters.config.getString( "abams.dataFolder" )
            + "AllenCatalogEnergies.matrix.cache";

    public static final String minerEnergyFilename = SetupParameters.config.getString( "abams.dataFolder" )
            + "AllenMinerEnergies.matrix.cache";
    public static final String minerMeanEnergyFilename = SetupParameters.config.getString( "abams.dataFolder" )
            + "AllenMinerMeanEnergies.matrix.cache";

    /**
     * @param args
     */
    public AllenCatalogMatrices() {
        this( false );
    }

    public AllenCatalogMatrices( boolean load ) {
        // if load is false then get it from the files
        if ( !load ) {
            try {
                readIn();
            } catch ( Exception e ) {
                e.printStackTrace();
                System.out.println( "Error loading matrices, loading from original csv files" );
                load = true;
            }
        }
        if ( load ) {
            try {
                loadFromOriginal();
                writeOut();
            } catch ( Exception e ) {
                throw new RuntimeException( e );
            }
        }

    }

    public void loadMinerOriginal() throws Exception {
        AllenMinerLoader miner = new AllenMinerLoader();
        log.info( "Loading allen miner" );
        minerEnergies = miner.getMatrix( AllenMinerLoader.totalExpressionInROI );

        minerMeanEnergies = miner.getMatrix( AllenMinerLoader.meanExpressionPosition );

    }

    public void loadFromOriginal() throws Exception {
        loadMinerOriginal();

        StructureExpressionInformationLoader loader = new StructureExpressionInformationLoader();
        log.info( "Loading energies" );
        energies = loader.getMatrix( StructureExpressionInformationLoader.energyPosition );

        log.info( "Loading levels" );
        levels = loader.getMatrix( StructureExpressionInformationLoader.levelPosition );

        log.info( "Loading densities" );
        densities = loader.getMatrix( StructureExpressionInformationLoader.densityPosition );

    }

    public void writeOut() throws IOException {
        ObjectOutputStream out = new ObjectOutputStream( new FileOutputStream( densityFilename ) );
        out.writeObject( densities );
        out.close();

        out = new ObjectOutputStream( new FileOutputStream( levelFilename ) );
        out.writeObject( levels );
        out.close();

        out = new ObjectOutputStream( new FileOutputStream( energyFilename ) );
        out.writeObject( energies );
        out.close();

        out = new ObjectOutputStream( new FileOutputStream( minerEnergyFilename ) );
        out.writeObject( minerEnergies );
        out.close();

        out = new ObjectOutputStream( new FileOutputStream( minerMeanEnergyFilename ) );
        out.writeObject( minerMeanEnergies );
        out.close();
    }

    public void readIn() throws Exception {
        ObjectInputStream in = new ObjectInputStream( new FileInputStream( densityFilename ) );
        densities = ( DenseDoubleMatrix<String, String> ) in.readObject();
        in.close();

        in = new ObjectInputStream( new FileInputStream( levelFilename ) );
        levels = ( DenseDoubleMatrix<String, String> ) in.readObject();
        in.close();

        in = new ObjectInputStream( new FileInputStream( energyFilename ) );
        energies = ( DenseDoubleMatrix<String, String> ) in.readObject();
        in.close();

        in = new ObjectInputStream( new FileInputStream( minerEnergyFilename ) );
        minerEnergies = ( DenseDoubleMatrix<String, String> ) in.readObject();
        in.close();

        in = new ObjectInputStream( new FileInputStream( minerMeanEnergyFilename ) );
        minerMeanEnergies = ( DenseDoubleMatrix<String, String> ) in.readObject();
        in.close();

        log.info( "Miner Energy Rows:" + minerEnergies.rows() );
        log.info( "Miner Energy Cols:" + minerEnergies.columns() );

        log.info( "Allen Energy Rows:" + energies.rows() );
        log.info( "Allen Energy Cols:" + energies.columns() );
    }

    public static void main( String[] args ) throws Exception {
        AllenCatalogMatrices matrices = new AllenCatalogMatrices( false );
        // matrices.loadMinerOriginal();
        // matrices.writeOut();

        // Util.writeRTable( "energies.R.table", matrices.getEnergies() );
    }

    public DoubleMatrix<String, String> getEnergies() {
        return energies;
    }

    public DoubleMatrix<String, String> getMinerEnergies() {
        return minerEnergies;
    }

    public DoubleMatrix<String, String> getMinerMeanEnergies() {
        return minerMeanEnergies;
    }

    public DoubleMatrix<String, String> getDensities() {
        return densities;
    }

    public DoubleMatrix<String, String> getLevels() {
        return levels;
    }

}
