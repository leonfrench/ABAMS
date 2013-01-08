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
package ubic.BAMSandAllen.gene2region;

import java.util.LinkedList;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.RegressionVector;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.math.CorrelationStats;
import ubic.basecode.math.Distance;
import cern.colt.list.DoubleArrayList;
import flanagan.analysis.Regression;

public class GeneToRegionCreator {
    private static Log log = LogFactory.getLog( GeneToRegionCreator.class.getName() );

    DoubleMatrix<String, String> dataMatrix;
    String gene;
    String region;

    // RegressionVector

    public GeneToRegionCreator( String geneName, String regionName, DoubleMatrix<String, String> expressionMatrix,
            DoubleMatrix<String, String> connectivityMatrix, DoubleMatrix<String, String> explainMatrix,
            RegressMatrix regressType ) throws Exception {
        this.gene = geneName;
        this.region = regionName;
        String regress = null;
        if ( regressType.equals( RegressMatrix.BOTH ) ) throw new RuntimeException( "Not supported" );
        if ( regressType.equals( RegressMatrix.EXPRESSION ) ) regress = gene;
        if ( regressType.equals( RegressMatrix.CONNECTIVITY ) ) regress = region;

        if ( !expressionMatrix.getColNames().equals( connectivityMatrix.getColNames() ) )
            throw new RuntimeException( "Matrix columns do not match" );
        if ( !expressionMatrix.getColNames().equals( explainMatrix.getColNames() ) )
            throw new RuntimeException( "Matrix columns do not match" );

        List<String> colNames = new LinkedList<String>( expressionMatrix.getColNames() );

        // leave out connections to itself
        if ( colNames.contains( region ) ) colNames.remove( region );

        dataMatrix = new DenseDoubleMatrix<String, String>( 3, colNames.size() );
        dataMatrix.setColumnNames( colNames );
        dataMatrix.addRowName( region );
        dataMatrix.addRowName( gene );
        dataMatrix.addRowName( "distance" );
        for ( String colName : colNames ) {
            dataMatrix.setByKeys( gene, colName, expressionMatrix.getByKeys( gene, colName ) );
            dataMatrix.setByKeys( region, colName, connectivityMatrix.getByKeys( region, colName ) );
            dataMatrix.setByKeys( "distance", colName, explainMatrix.getByKeys( region, colName ) );
        }

        Regression regression = new Regression( dataMatrix.getRowByName( "distance" ), dataMatrix
                .getRowByName( regress ) );
        // regression.supressErrorMessages();
        // regression.supressPrint();
        // regression.supressStats();
        if ( Util.countValues( dataMatrix, Double.NaN ) > 0 ) {
            log.info( "Skipping regression" );

        } else {
            //log.info( "Doing regression" );
            regression.linear();
            // replace target with residuals
            double[] residuals = regression.getResiduals();
            int regressIndex = dataMatrix.getRowIndexByName( regress );
            for ( int i = 0; i < colNames.size(); i++ ) {
                dataMatrix.set( regressIndex, i, residuals[i] );
            }
            //Util.writeRTable( "/grp/java/workspace/BAMSandAllen/data/incomplete.txt", dataMatrix );
        }

    }

    public GeneToRegionCreator( String geneName, String regionName, DoubleMatrix<String, String> expressionMatrix,
            DoubleMatrix<String, String> connectivityMatrix ) {

        this.gene = geneName;
        this.region = regionName;

        if ( !expressionMatrix.getColNames().equals( connectivityMatrix.getColNames() ) )
            throw new RuntimeException( "Matrix columns do not match" );

        List<String> colNames = new LinkedList<String>( expressionMatrix.getColNames() );

        // leave out connections to itself
        if ( colNames.contains( regionName ) ) colNames.remove( regionName );

        dataMatrix = new DenseDoubleMatrix<String, String>( 2, colNames.size() );
        dataMatrix.setColumnNames( colNames );
        dataMatrix.addRowName( region );
        dataMatrix.addRowName( gene );
        for ( String colName : colNames ) {
            dataMatrix.setByKeys( gene, colName, expressionMatrix.getByKeys( geneName, colName ) );
            dataMatrix.setByKeys( region, colName, connectivityMatrix.getByKeys( regionName, colName ) );
        }
    }

    public double getCorrelation() {
        return CorrelationStats.correl( dataMatrix.getRowByName( gene ), dataMatrix.getRowByName( region ) );
    }

    public double getJackKnifeCorrelation() {
        double totalCorrelation = 0;
        for ( int i = 0; i < dataMatrix.columns(); i++ ) {
            DoubleArrayList exp = new DoubleArrayList( dataMatrix.getRowByName( gene ) );
            DoubleArrayList con = new DoubleArrayList( dataMatrix.getRowByName( region ) );
            exp.remove( i );
            con.remove( i );
            totalCorrelation += CorrelationStats.correl( exp.elements(), con.elements() );
        }
        totalCorrelation = totalCorrelation / dataMatrix.columns();
        return totalCorrelation;
    }

    public double getRankCorrelation() {
        return Distance.spearmanRankCorrelation( new DoubleArrayList( dataMatrix.getRowByName( gene ) ),
                new DoubleArrayList( dataMatrix.getRowByName( region ) ) );
    }

    public String toString() {
        return dataMatrix.toString();
    }
}
