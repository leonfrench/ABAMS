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
package ubic.BAMSandAllen.ClassSelectors;

import ubic.BAMSandAllen.Vocabulary;

import com.hp.hpl.jena.ontology.OntClass;

public class Major17ClassSelector extends BrainRegionClassSelector implements ClassSelector {
    public boolean test( OntClass ontClass ) {
        if (!super.test(ontClass)) {
            return false;
        }
        return ontClass.hasProperty( Vocabulary.hasAllen17Mapping );
    }
}
