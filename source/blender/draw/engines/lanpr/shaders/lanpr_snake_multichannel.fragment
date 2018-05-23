//#version 130
smooth in vec3 normal;
smooth in vec4 finalcolor;

float Interpolate(float between1,float between2,float value1,float value2,float key){
    float i = (key-between1)/(between2-between1);
    return value1*(1-i)+value2*i;
}

void main(){
    float value = dot(vec3(0,0,1),normal);
    //if(value<0.65) value=0.15;
    //else if(value>=0.65 && value<0.85) value=Interpolate(0.65,0.85,0.15,0.75,value);
    //else if(value>=0.85 && value<0.95) value=0.75;
    //else if(value>=0.95) value=0.9;
    gl_FragData[0] = vec4(vec3(finalcolor)*value,1);
    gl_FragData[1] = vec4(normal,1);//vec4((normal+vec3(1))*0.5,1);
}