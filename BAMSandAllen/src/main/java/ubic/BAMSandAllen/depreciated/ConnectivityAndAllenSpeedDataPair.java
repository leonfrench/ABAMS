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
package ubic.BAMSandAllen.depreciated;

import java.util.LinkedList;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.LOOCorrelObject;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenExpressionMatrixPair;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.math.CorrelationStats;

/**
 * This class allows leave one out correlation computation for the matrix B
 */
@Deprecated
public class ConnectivityAndAllenSpeedDataPair {// extends ConnectivityAndAllenExpressionMatrixPair {
// private static Log log = LogFactory.getLog( ConnectivityAndAllenSpeedDataPair.class.getName() );
// // correlation object corresponding to the triangles
// LOOCorrelObject[] corrTriangleA;
// LOOCorrelObject[] corrTriangleB;
//
// public ConnectivityAndAllenSpeedDataPair( BrainRegionClassSelector selector, boolean doLog, boolean onePlusLog,
// boolean square, double zeroReplacement, String matrixName, Direction direction ) throws Exception {
// super( selector, doLog, onePlusLog, square, zeroReplacement, matrixName, direction );
//
// }
//
// public double[] getLOOCorrelationTriangleB( String removed ) {
// // add it back, don't remove it permenantly
// return getLOOCorrelationTriangleB( removed, true );
// }
//
// public double[] getLOOCorrelationTriangleB( String removed, boolean addBack ) {
// double[] result = new double[corrTriangleB.length];
// double[] xj = new double[corrTriangleB.length];
// double[] yj = new double[corrTriangleB.length];
//
// int row = matrixB.getRowIndexByName( removed );
// int size = matrixB.columns();
//
// int spot = 0;
// for ( int i = 0; i < size; i++ ) {
// for ( int j = i + 1; j < size; j++ ) {
// xj[spot] = matrixB.get( row, i );
// yj[spot] = matrixB.get( row, j );
// spot++;
// }
// }
//
// for ( int i = 0; i < corrTriangleB.length; i++ ) {
// corrTriangleB[i].removePair( xj[i], yj[i] );
// result[i] = corrTriangleB[i].correl();
// // put it back
// if ( addBack ) {
// corrTriangleB[i].addPair( xj[i], yj[i] );
// }
// }
//
// return result;
// }
//
// public double[] getCorrelationTriangle( DoubleMatrix<String, String> matrix ) {
// if ( matrix.equals( matrixB ) ) {
// // use the correlation triangle for speed
// double[] result = new double[corrTriangleB.length];
// for ( int i = 0; i < result.length; i++ )
// result[i] = corrTriangleB[i].correl();
// return result;
// }
// if ( matrix.equals( matrixA ) ) {
// // use the correlation triangle for speed
// double[] result = new double[corrTriangleA.length];
// for ( int i = 0; i < result.length; i++ )
// result[i] = corrTriangleA[i].correl();
// return result;
// }
// // should be a or b
// return null;
// }
//
// public double correWithoutRow( String rowName ) {
// double real = CorrelationStats
// .correl( getCorrelationTriangle( matrixA ), getLOOCorrelationTriangleB( rowName ) );
// return real;
// }
//
// public LOOCorrelObject[] setupCorrelationTriangle( DoubleMatrix<String, String> matrix ) {
// int size = matrix.columns();
// LOOCorrelObject[] resultTriangle = new LOOCorrelObject[size * ( size - 1 ) / 2];
//
// int spot = 0;
// for ( int i = 0; i < size; i++ ) {
// for ( int j = i + 1; j < size; j++ ) {
// resultTriangle[spot] = new LOOCorrelObject();
// resultTriangle[spot].correlAll( matrix.getColumn( i ), matrix.getColumn( j ) );
// spot++;
// }
// }
// return resultTriangle;
// }
//
// public void setupCorrelationTriangles() {
// corrTriangleA = setupCorrelationTriangle( matrixA );
// corrTriangleB = setupCorrelationTriangle( matrixB );
// }
//
// // make sure we update the triangle when columns are changed
// public void sameSpace() {
// super.sameSpace();
// setupCorrelationTriangles();
// }
//
// public void slimMatricesOnce() {
// super.slimMatricesOnce();
// setupCorrelationTriangles();
// }
//
// public void removeConnectionZeroes() {
// super.removeConnectionZeroes();
// setupCorrelationTriangles();
// }
//
// public void removeBedNucleiStria() {
// super.removeBedNucleiStria();
// setupCorrelationTriangles();
// }
//
// public void removeDataRow( String row ) {
// super.removeDataRow( row );
// // this calls removeDataRows
// // setupCorrelationTriangles();
// }
//
// // remove it fast, two layers, edit the matrix and the triangle
// public void removeDataRows( List<String> rows ) {
// // remove from the triangle correlations first
// for ( String row : rows ) {
// boolean addBack = false;
// getLOOCorrelationTriangleB( row, addBack );
// }
//
// // remove from the matrix
// super.removeDataRows( rows );
// // setupCorrelationTriangles();
// }
//
// public void squareConnectionMatrix() {
// super.squareConnectionMatrix();
// setupCorrelationTriangles();
// }
//
// // setup corrObject
// // get triangle size
// // make corrObject for each spot
//
// // get corr
//
// // remove row and get corr - remove it from every spot in triangle, and underlying matrix? optional
//
// public static void main( String[] args ) throws Exception {
// // below is testing code
// long t;
// ConnectivityAndAllenExpressionMatrixPair forR = new ConnectivityAndAllenExpressionMatrixPair(
// new BrainRegionClassSelector(), true, false, false, Double.NaN, "NewEnergies",
// AnalyzeBAMSandAllenGenes.Direction.OUTGOING );
// forR.run();
// log.info( forR.getCorrelation() );
//
// ConnectivityAndAllenSpeedDataPair forR2 = new ConnectivityAndAllenSpeedDataPair(
// new BrainRegionClassSelector(), true, false, false, Double.NaN, "NewEnergies",
// AnalyzeBAMSandAllenGenes.Direction.OUTGOING );
// forR2.run();
// // forR2.getCorrelationTriangle( null );
// t = System.currentTimeMillis();
// log.info( forR.getCorrelation( true ) );
// log.info( forR.getCorrelation( true ) );
// log.info( forR.getCorrelation( true ) );
// log.info( System.currentTimeMillis() - t );
//
// t = System.currentTimeMillis();
// log.info( forR2.getCorrelation( true ) );
// log.info( forR2.getCorrelation( true ) );
// log.info( forR2.getCorrelation( true ) );
// log.info( System.currentTimeMillis() - t );
//
// String row = forR2.getMatrixBDataRows().get( 4 );
//
// log.info( row );
// List<String> rows = new LinkedList<String>();
// rows.add( row );
//
// t = System.currentTimeMillis();
// log.info( "REFACTORED:" );
// log.info( forR.correWithoutMatrixBDataRow( row ) );
// log.info( forR.correWithoutMatrixBDataRow( row ) );
// log.info( forR.correWithoutMatrixBDataRow( row ) );
//
// log.info( System.currentTimeMillis() - t );
//
// t = System.currentTimeMillis();
// log.info( "SLOW:" );
// // log.info( forR.correlationReducedDataMatrix( rows ) );
// // log.info( forR.correlationReducedDataMatrix( rows ) );
// log.info( forR.correlationReducedDataMatrix( rows ) );
// log.info( forR.correlationReducedDataMatrix( rows ) );
// log.info( forR.correlationReducedDataMatrix( rows ) );
//
// t = System.currentTimeMillis();
// log.info( "FAST:" );
// log.info( forR2.correWithoutRow( row ) );
// log.info( forR2.correWithoutRow( row ) );
// log.info( forR2.correWithoutRow( row ) );
// log.info( System.currentTimeMillis() - t );
//
// log.info( forR2.getCorrelation( true ) );
//
// }
//
}
