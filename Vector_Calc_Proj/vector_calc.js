// GET JS MATRIX DATA
var a = JSMatrix.Vector.create([1, 2, 3]);

/******************* PUT EVERYTHING IN ONE BIG FUNCTION AND USE EVENT LISTENERS TO CALL SUBROUTINES *********************/

/*
// ** Gathering Inputs: create inputs and store values in arrays to pull from

var angles = []; // vector sized by # of members, each element is a member's angle
var text_vector = [];
var member_names = [];
var number = 0; // number of members
var ff = 1; // incremnter (why is it 1?)

var joint_text_vector = [];
var joints = []; 
var jj = 1; // incremnter (why is it 1?)

var appliedForces = [];
var aF = 1; // incremnter (why is it 1?)
var hmfExecuted = false;
// forces will increase drop down options for joints, joints will use drop downs to select forces they
*/

function execute(){
    var angles = []; // vector sized by # of members, each element is a member's angle
    var member_names = [];
    var number = 0; // number of members

    var joints = [];
    var number_j = 0; 

    // var appliedForces = [];

    howManyForces(number, member_names, angles, joints, number_j); // functions are nested to make sure variables mutate
}

function howManyForces(number, member_names, angles, joints, number_j) {
    // Number of inputs to create
    number = document.getElementById("forces").value;
    if (number == 0 || number == ''){
        alert('Input a number of Members/Reactions greater than 0');
        return;
    }
    // Container <div> where dynamic content will be placed
    var container = document.getElementById("memberFormContainer");
    // Clear previous contents of the container
    while (container.hasChildNodes()) { // make sure not to remove the original child: use if-statement
        container.removeChild(container.lastChild);
    }
    for (i=0;i<number;i++){
        // Append a node with a random text
        container.appendChild(document.createTextNode("Member " + (i+1)));
        // Create an <input> element, set its type and name attributes
        var input = document.createElement("input");
        input.type = "number";
        input.step = 0.00000000001;
        input.id = "member" + i;
        member_names.push("member" + i); //save member names in member_names vector
        container.appendChild(input);
        // Append a line break 
        container.appendChild(document.createElement("br"));
    }
    document.getElementById('memberButton2').addEventListener('click', function() { submitForces(number, member_names, angles, joints, number_j);
    }, false); // watch what happens to the values called in the event listener
}

/** make a function called submitForces() that responds to the submit button under the members section:
 * if howManyForces() hasn't been executed --> error message, exit
    * while container isn't empty:
    *      if a field is empty --> error message, exit
    *      if a field has negative value --> error message, exit
    *      if a field has value greater than 180 --> error message, exit
    *      else: 
    *           take field value and push it in the back of forces vector
 * **/
function submitForces(number, member_names, angles, joints, number_j){
    angles = [];
    document.getElementById('memberButton2').style.backgroundColor = "blue";
    // make sure function isn't called too early   
    for (i=0;i<number;i++){ 
        var value = document.getElementById(member_names[i]).value;
        if (value == "") { // checked!
            document.getElementById('memberButton2').style.backgroundColor = "gray";
            angles = [];
            alert('One of your fields are empty'); // doesnt show up
            return; //break instead?
        }
        if ( (isNaN(value)) || (value < 0) || (value > 180) ){ // checked!
            document.getElementById('memberButton2').style.backgroundColor = "gray";
            angles = [];
            alert('Please Input a positive number between 0-180'); // doesnt show up
            return;
        }
        angles.push(value); // save angle value of member in 'angles' vector
    }
    document.getElementById('jointButton').addEventListener('click', function() { howManyJoints(number_j, number, joints, angles);
    }, false);
} 

function howManyJoints(number_j, number, joints, angles){
    // DO NOT EXECUTE UNLESS MEMBERS HAVE BEEN SUBMITTED
    // Number of inputs to create
    number_j = document.getElementById("joints").value;
    // Container <div> where dynamic content will be placed
    var container = document.getElementById("jointFormContainer");
    // Clear previous contents of the container
    while (container.hasChildNodes()) { 
        container.removeChild(container.lastChild);
    }
    for (i=0;i<number_j;i++){
        // joint_text_vector.push("Joint " + (i+1)); //save text_node value in text_vector
        // Append a node with a random text
        container.appendChild(document.createTextNode("Joint " + (i+1)));
        // Create an array of checkbox <input> elements, set type, id attributes
        for (j=0;j<number;j++){
            var input = document.createElement("input");
            input.type = "checkbox";
            input.id = "joint" + i + "member" + j;
            container.appendChild(input);
        }
        joints.push("joint" + i); //save joint names in joints vector
        // Append a line break 
        container.appendChild(document.createElement("br"));
    }
    document.getElementById('jointButton2').addEventListener('click', function() { submitJoints(number_j, number, angles);
    }, false);
}

/**
 * Create joint "connections" matrix --> [joints, members]
 * have a 1 for each checked member of a joint and a 0 for each non checked
 */

function submitJoints(number_j, number, angles){
    document.getElementById('jointButton2').style.backgroundColor = "blue";
    var connection = new JSMatrix.Matrix(number_j, number);
    for (i=0;i<number_j;i++){
        for (j=0;j<number;j++){
           var joint_member = document.getElementById("joint" + i + "member" + j);
            if (joint_member.checked){
                connection.set(i,j,1);
            }
        }
    }
    document.getElementById("appliedForceButton").addEventListener('click', function() { howManyaddAppliedForces(connection, number_j, angles, number);
    }, false);
}

/**
 * 
 */
function howManyaddAppliedForces(connection, number_j, angles, number){
    // Number of inputs to create
   // number = document.getElementById("forces").value;
    // Container <div> where dynamic content will be placed
    var container = document.getElementById("appliedForceContainer");
    // Clear previous contents of the container
    while (container.hasChildNodes()) { // make sure not to remove the original child: use if-statement
        container.removeChild(container.lastChild);
    }
    for (i=0;i<number_j;i++){
        // Append a node with a random text
        container.appendChild(document.createTextNode("Joint " + (i+1) + " " + "Fx"));
        // Create an <input> element, set its type and name attributes
        var input = document.createElement("input");
        input.type = "number";
        input.step = 0.00000000001;
        input.id = "aFx" + i;
        // member_names.push("member" + i); //save member names in member_names vector
        container.appendChild(input);
        // Append a line break 
        container.appendChild(document.createElement("br"));

        container.appendChild(document.createTextNode("Joint " + (i+1) + " " + "Fy"));
        var input2 = document.createElement("input");
        input2.type = "number";
        input2.step = 0.00000000001;
        input2.id = "aFy" + i;
        container.appendChild(input2);
        container.appendChild(document.createElement("br"));
    }
    document.getElementById('calculate').addEventListener('click', function() { eqnSolver(connection, angles, number_j, number);
    }, false); // watch what happens to the values called in the event listener
}


// Building Calculator

function dotProduct(connectionRow, angleRow){
    var newRow = new JSMatrix.Vector(connectionRow.size());
    for (e = 0; e<connectionRow.size(); e++){
        newRow.set(e, (connectionRow.get(e) * angleRow.get(e)));
    }
    return newRow;
}

function cosD(anglesVector){ // this function changes variables it shouldn't be
    for (e = 0; e<anglesVector.size(); e++){ // convert degrees to radians
        anglesVector.set(e, (anglesVector.get(e) * (Math.PI/180)));
        anglesVector.set(e, Math.cos(anglesVector.get(e))); //get cos value for each radian
    }
    return anglesVector;
}

function sinD(anglesVector){
    for (e = 0; e<anglesVector.size(); e++){ // convert degrees to radians
        anglesVector.set(e, (anglesVector.get(e) * (Math.PI/180)));
        anglesVector.set(e, Math.sin(anglesVector.get(e))); //get cos value for each radian
    }
    return anglesVector;
}

function getAppliedForces(number_j){
    var aFVector = [];
    for(i=0; i<number_j; i++){
        aFVector.push(document.getElementById('aFx' + i).value);
        aFVector.push(document.getElementById('aFy' + i).value);
    }
    return aFVector;
}

function eqnSolver(connection, angles, number_j, number){
    document.getElementById('calculate').style.backgroundColor = "blue";
    for(i=0; i<number_j; i++){
        if (document.getElementById('aFx' + i).value == '' || document.getElementById('aFy' + i).value == ''){
            aFVector = [];
            alert('One of your Applied Force Fields are Empty');
            return;
        }
    }
    var A = new JSMatrix.Matrix(number_j*2, number);
    var angles_vector = JSMatrix.Vector.create(angles); // may need to use angles to create matrix or vector
    var aF = getAppliedForces(number_j);
    var aF_Vector = JSMatrix.Vector.create(aF);
    for (i = 0; i<number_j; i++){
        for (j = 0; j<2; j++){ 
            var row_num = (2*i)+j;
            if (j == 0) {
                // copy "A" actions here
                var new_angles_vector = new JSMatrix.Vector(angles_vector.size());
                for (e = 0; e<angles_vector.size(); e++){
                    new_angles_vector.set(e, angles_vector.get(e));
                }
                var newRow = dotProduct(connection.getRow(i), cosD(new_angles_vector));
                A.setRow(row_num, newRow);
            }
            else {
                // copy "A" actions here
                var new_angles_vector_b = new JSMatrix.Vector(angles_vector.size());
                for (e = 0; e<angles_vector.size(); e++){
                    new_angles_vector_b.set(e, angles_vector.get(e));
                }
                var newRow = dotProduct(connection.getRow(i), sinD(new_angles_vector_b));
                A.setRow(row_num, newRow);
            }
        }
        angles_vector.iadd(connection.getRow(i).mul(180)); // negate the angle vector's radians but only the ones that were in the previous connection
    }
    var answers = A.linSolve(aF_Vector.mul(-1));
    
    printAnswer(answers);
}

function printAnswer(answersVector){
    var container = document.getElementById('answersContainer');
    for (i = 0; i < answersVector.size(); i++){
        container.appendChild(document.createTextNode("Member " + (i+1) + ": " + answersVector[i] + " units of Force"));
        container.appendChild(document.createElement("br"));
    }
}

function inputDefaultValues(){
    document.getElementById('defaultValues').style.backgroundColor = "blue";
    document.getElementById('memberButton2').style.backgroundColor = "blue";
    document.getElementById('jointButton2').style.backgroundColor = "blue";
    document.getElementById('memberButton').style.backgroundColor = "blue";
    document.getElementById('jointButton').style.backgroundColor = "blue";
    document.getElementById('appliedForceButton').style.backgroundColor = "blue";

    var container1 = document.getElementById('memberFormContainer');
    var container2 = document.getElementById('jointFormContainer');
    var container3 = document.getElementById('appliedForceContainer');
    var memberAngles = [40, 62.62, 130, 62.28, 40, 90, 0, 0];
    var connection = [ [1, 1, 0, 0, 0, 0, 0, 0], [1, 0, 1, 1, 0, 0, 0, 1], [0, 1, 1, 0, 1, 0, 0, 0], [0, 0, 0, 1, 1, 1, 1, 0] ];
    var connectionMatrix = JSMatrix.Matrix.create(connection);
    var aF = [0, -8829, 0, 0, 0, 0, 0, 0];


    if (container3.childElementCount > 1){
        while (container3.hasChildNodes()) { // make sure not to remove the original child: use if-statement
            container3.removeChild(container3.lastChild);
        }
    }


    if (container2.childElementCount > 2){
        while (container2.hasChildNodes()) { // make sure not to remove the original child: use if-statement
            container2.removeChild(container2.lastChild);
        }
    }
    
    if (container1.childElementCount > 2){
        while (container1.hasChildNodes()) { // make sure not to remove the original child: use if-statement
            container1.removeChild(container1.lastChild);
        }
    }

    // fill Members
    for (i=0;i<8;i++){
        // Append a node with a random text
        container1.appendChild(document.createTextNode("Member " + (i+1)));
        // Create an <input> element, set its type and name attributes
        var input1 = document.createElement("input");
        input1.type = "number";
        input1.step = 0.00000000001;
        input1.id = "member" + i;
        input1.value = memberAngles[i];
        container1.appendChild(input1);
        // Append a line break 
        container1.appendChild(document.createElement("br"));
    }

    for (i=0;i<4;i++){
        // joint_text_vector.push("Joint " + (i+1)); //save text_node value in text_vector
        // Append a node with a random text
        container2.appendChild(document.createTextNode("Joint " + (i+1)));
        // Create an array of checkbox <input> elements, set type, id attributes
        for (j=0;j<8;j++){
            var input2 = document.createElement("input");
            input2.type = "checkbox";
            input2.id = "joint" + i + "member" + j;
            if (connection[i][j] == 1){
                input2.checked = true;
            }
            container2.appendChild(input2);
        }
        // Append a line break 
        container2.appendChild(document.createElement("br"));
    }

    for (i=0;i<4;i++){
        // Append a node with a random text
        container3.appendChild(document.createTextNode("Joint " + (i+1) + " " + "Fx"));
        // Create an <input> element, set its type and name attributes
        var input3 = document.createElement("input");
        input3.type = "number";
        input3.step = 0.00000000001;
        input3.id = "aFx" + i;
        // member_names.push("member" + i); //save member names in member_names vector
        container3.appendChild(input3);
        // Append a line break 
        container3.appendChild(document.createElement("br"));

        container3.appendChild(document.createTextNode("Joint " + (i+1) + " " + "Fy"));
        var input4 = document.createElement("input");
        input4.type = "number";
        input4.step = 0.00000000001;
        input4.id = "aFy" + i;
        container3.appendChild(input4);
        container3.appendChild(document.createElement("br"));
    }

    for (i=0;i<4;i++){ // may have to refer backwards instead of forward
        var jj = i*2;
        var kk = (i*2)+1;
        document.getElementById("aFy" + i).value = aF[kk];
        document.getElementById("aFx" + i).value = aF[jj];
    }


    document.getElementById('calculate').addEventListener('click', function() { eqnSolver(connectionMatrix, memberAngles, 4, 8);
    }, false);
}



