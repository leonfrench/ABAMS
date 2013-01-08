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

import com.hp.hpl.jena.ontology.OntClass;
import com.hp.hpl.jena.rdf.model.Property;
import com.hp.hpl.jena.rdf.model.RDFNode;
import com.hp.hpl.jena.rdf.model.Resource;
import com.hp.hpl.jena.rdf.model.Selector;
import com.hp.hpl.jena.rdf.model.Statement;

public class CombinedSelector implements Selector {
    Selector[] selectors;

    public CombinedSelector(Selector... selectors  ) {
        this.selectors = selectors;
    }

    public RDFNode getObject() {
        return null;
    }

    public Property getPredicate() {
        // TODO Auto-generated method stub
        return null;
    }

    public Resource getSubject() {
        // TODO Auto-generated method stub
        return null;
    }

    public boolean isSimple() {
        boolean result = true;
        for (Selector s : selectors) {
            result = result && s.isSimple();
        }
        return result;
    }

    public boolean test( Statement stmt ) {
        boolean result = true;
        for (Selector s : selectors) {
            result = result && s.test( stmt );
        }
        return result;
    }

}
