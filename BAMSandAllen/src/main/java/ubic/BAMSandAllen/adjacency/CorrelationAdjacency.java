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
package ubic.BAMSandAllen.adjacency;

import java.io.Serializable;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.LOOCorrelObject;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.MatrixPairs.MatrixPair;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class CorrelationAdjacency implements AdjacencyCompute, Serializable {
    private static Log log = LogFactory.getLog( CorrelationAdjacency.class.getName() );

    boolean spearman;
    LOOCorrelObject[][] corrTriangle;
    DoubleMatrix<String, String> matrix;
    boolean setup;
    DoubleMatrix<String, String> resultPreConstructed;

    DoubleMatrix<String, String> fullAdjacency;

    public CorrelationAdjacency( DoubleMatrix<String, String> matrix ) {
        this( matrix, false );
        setup = false;

    }

    public CorrelationAdjacency( DoubleMatrix<String, String> matrixP, boolean spearman ) {
        this.spearman = spearman;
        setMatrix( matrixP, true );

    }

    public CorrelationAdjacency clone() {
        CorrelationAdjacency result = new CorrelationAdjacency( matrix, spearman );
//        result.setup = setup;
//        int size = matrix.columns();
//        result.corrTriangle = new LOOCorrelObject[size][];
//
//        if ( corrTriangle != null ) {
//            for ( int i = 0; i < size; i++ ) {
//                result.corrTriangle[i] = new LOOCorrelObject[i];
//                for ( int j = 0; j < i; j++ ) {
//                    result.corrTriangle[i][j] = corrTriangle[i][j].clone();
//                }
//            }
//        }
//
//        result.resultPreConstructed = resultPreConstructed;
//        result.fullAdjacency = fullAdjacency;

        return result;
    }

    public void setMatrix( DoubleMatrix<String, String> matrix, boolean safe ) {
        // log.info( "Calling setmatrix size:" + matrix.rows() + " x " + matrix.columns() );
        this.matrix = matrix;

        // init these matrices once instead of everytime
        resultPreConstructed = new DenseDoubleMatrix<String, String>( matrix.columns(), matrix.columns() );
        resultPreConstructed.setRowNames( matrix.getColNames() );
        resultPreConstructed.setColumnNames( matrix.getColNames() );

        if ( safe ) setupAdjacency();
    }

    /*
     * assumes same space (non-Javadoc)
     * 
     * @see ubic.BAMSandAllen.refactored.TriangleCompute#setupCorrelationTriangle()
     */
    public void setupAdjacency() {
        setup = false;
    }

    private void actuallySetupCorrelationTriangle() {
        log.info( "Setting up triangle" );
        int size = matrix.columns();
        corrTriangle = new LOOCorrelObject[size][];

        // fill the bottom right triangle
        for ( int i = 0; i < size; i++ ) {
            corrTriangle[i] = new LOOCorrelObject[i];
            for ( int j = 0; j < i; j++ ) {
                corrTriangle[i][j] = new LOOCorrelObject();
                corrTriangle[i][j].correlAll( matrix.getColumn( i ), matrix.getColumn( j ) );
            }
        }
        setup = true;
    }

    public DoubleMatrix<String, String> getAdjacency() {

        if ( !setup ) actuallySetupCorrelationTriangle();
        // return the cached version if we can
        if ( fullAdjacency == null ) {
            fullAdjacency = new DenseDoubleMatrix<String, String>( matrix.columns(), matrix.columns() );
            fullAdjacency.setRowNames( matrix.getColNames() );
            fullAdjacency.setColumnNames( matrix.getColNames() );

            for ( int i = 0; i < corrTriangle.length; i++ ) {
                for ( int j = 0; j < i; j++ ) {
                    fullAdjacency.set( i, j, corrTriangle[i][j].correl() );
                    fullAdjacency.set( j, i, corrTriangle[i][j].correl() );
                }
            }
        }
        return fullAdjacency;
    }

    // public double[] getTriangle() {
    // if ( !setup ) actuallySetupCorrelationTriangle();
    //
    // int size = corrTriangle.length;
    // double[] result = new double[size * ( size - 1 ) / 2];
    // int spot = 0;
    // for ( int i = 0; i < size; i++ ) {
    // for ( int j = 0; i < j; j++ ) {
    // result[spot++] = corrTriangle[i][j].correl();
    // }
    // }
    // return result;
    // }

    public DoubleMatrix<String, String> getAdjacency( String removed ) {
        // add it back, don't remove it permenantly
        return getAdjacency( removed, true );
    }

    public DoubleMatrix<String, String> getAdjacency( String removed, boolean adback ) {
        if ( !setup ) actuallySetupCorrelationTriangle();
        // log.info( resultPreConstructed.columns() + "X" + resultPreConstructed.rows() );
        // if we are not putting it back we gota setup again (afterwards)

        int row = matrix.getRowIndexByName( removed );

        double xj;
        double yj;

        for ( int i = 0; i < corrTriangle.length; i++ ) {
            xj = matrix.get( row, i );
            for ( int j = 0; j < i; j++ ) {
                yj = matrix.get( row, j );

                corrTriangle[i][j].removePair( xj, yj );
                double value = corrTriangle[i][j].correl();
                // result.set( i, j, value ); // removed for speed
                resultPreConstructed.set( j, i, value );

            }
        }

        // put it back into the correlation triangle
        if ( adback ) {
            for ( int i = 0; i < corrTriangle.length; i++ ) {
                xj = matrix.get( row, i );
                for ( int j = 0; j < i; j++ ) {
                    yj = matrix.get( row, j );
                    corrTriangle[i][j].addPair( xj, yj );
                }
            }
        } else {
            fullAdjacency = null;
        }

        return resultPreConstructed;
    }
}
