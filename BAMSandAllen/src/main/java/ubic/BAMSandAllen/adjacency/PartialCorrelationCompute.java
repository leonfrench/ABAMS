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

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.RegressionVector;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import flanagan.analysis.Regression;

/**
 * Computes adjacency against the residuals of a regression.
 * 
 * @author leon
 * @version $Id
 */
public class PartialCorrelationCompute implements AdjacencyCompute, Serializable {
    private static Log log = LogFactory.getLog( PartialCorrelationCompute.class.getName() );
    final String targetString = "Target";
    final String explainString = "Explain";
    final String residualsString = "Residuals";
    final String calcString = "Calculated";

    ABAMSDataMatrix explainMatrix;
    DoubleMatrix<String, String> targetMatrix;
    DoubleMatrix<String, String> fullAdjacency;

    AdjacencyCompute computeBeforeResiduals;
    RegressionVector triangles;
    boolean computeRegression;

    public PartialCorrelationCompute( DoubleMatrix<String, String> targetMatrix, ABAMSDataMatrix subtract,
            AdjacencyCompute computeBeforeResiduals ) {

        this.explainMatrix = subtract;
        this.targetMatrix = targetMatrix;
        this.computeBeforeResiduals = computeBeforeResiduals;
        boolean safe = true;
        computeBeforeResiduals.setMatrix( targetMatrix, safe );

        triangles = new RegressionVector( 4, targetMatrix.columns() );
        triangles.add( explainString, explainMatrix.getAdjacency() );
        computeRegression = true;

        fullAdjacency = null;
        //log.info( "new matrix made:" + targetMatrix.rows() );
    }

    public ABAMSDataMatrix getExplainMatrix() {
        return explainMatrix;
    }

    public AdjacencyCompute clone() {
        // make sure the clone is passed on
        PartialCorrelationCompute result = new PartialCorrelationCompute( targetMatrix, explainMatrix,
                computeBeforeResiduals.clone() );
        return result;
    }

    public DoubleMatrix<String, String> getAdjacency() {
        // return the cache if we can
        if ( fullAdjacency == null ) {
            DoubleMatrix<String, String> targetAdjacency = computeBeforeResiduals.getAdjacency();
            fullAdjacency = getAdjacency( targetAdjacency );
        }
        return fullAdjacency;
    }

    public void setComputeRegression( boolean computeRegression ) {
        // compute it one more time before possibly turning it off (in case any were removed)
        if ( this.computeRegression && !computeRegression ) getAdjacency();
        this.computeRegression = computeRegression;
    }

    public boolean getComputeRegression(  ) {
        return computeRegression;
    }

    
    public void setTriangles( RegressionVector triangles ) {
        this.triangles = triangles;
    }

    public RegressionVector getTriangles() {
        return triangles;
    }

    public void reComputeCalculatedValues() {
        Regression regression = new Regression( triangles.getRow( explainString ), triangles.getRow( targetString ) );
        regression.linear();

        // triangles.add( residualsString, regression.getResiduals() );

        triangles.add( calcString, regression.getYcalc() );
    }

    public DoubleMatrix<String, String> getAdjacency( DoubleMatrix<String, String> targetAdjacency ) {

        // put them in same space
        triangles.add( targetString, targetAdjacency );

        double[] calculated;
        if ( computeRegression ) {
            reComputeCalculatedValues();
        }

        calculated = triangles.getRow( calcString );

        double[] target = triangles.getRow( targetString );

        int length = calculated.length;

        double residuals[] = new double[length];
        for ( int i = 0; i < length; i++ ) {
            residuals[i] = target[i] - calculated[i];
        }

        triangles.add( residualsString, residuals );

        DoubleMatrix<String, String> result = triangles.getRowAsMatrix( residualsString );

        return result;
    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {

    }

    // needed?
    public void setupAdjacency() {
        computeBeforeResiduals.setupAdjacency();
    }

    public void setMatrix( DoubleMatrix<String, String> matrix, boolean safe ) {
        targetMatrix = matrix;
        fullAdjacency = null;
        computeBeforeResiduals.setMatrix( matrix, safe );
    }

    public DoubleMatrix<String, String> getAdjacency( String removed, boolean adback ) {
        DoubleMatrix<String, String> targetAdjacency = computeBeforeResiduals.getAdjacency( removed, adback );
        if ( !adback ) {
            // make it recompute the adjacency
            fullAdjacency = null;
        }
        // / below will minus residuals
        return getAdjacency( targetAdjacency );
    }

    public DoubleMatrix<String, String> getAdjacency( String removed ) {
        return getAdjacency( removed, true );
    }
}
