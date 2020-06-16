from OpenGL.GL import *
import scipy.spatial
from klampt.vis import glcommon
from klampt.math import vectorops

def draw_hull(hull,scale=1.0):
    """Raw GL draw of scipy.spatial.ConvexHull object as a polytope.  Can be lit and/or 
    rendered with a 3D texture.

    The object is scaled by scale.
    """
    if len(hull.simplices)==0:
        print "Hull with no simplices?"
        return
    centroid = sum([hull.points[v] for v in hull.vertices],[0.0]*3) / len(hull.vertices) * scale
    vmin = [min([hull.points[v][i] for v in hull.vertices])*scale for i in range(3)]
    vmax = [max([hull.points[v][i] for v in hull.vertices])*scale for i in range(3)]
    uvwscale = [1.0/(b-a) if b != a else 1.0 for (a,b) in zip(vmin,vmax)]

    for simplex in hull.simplices:
        ia,ib,ic = simplex
        a = vectorops.mul(hull.points[ia].tolist(),scale)
        b = vectorops.mul(hull.points[ib].tolist(),scale)
        c = vectorops.mul(hull.points[ic].tolist(),scale)
        n = vectorops.cross(vectorops.sub(b,a),vectorops.sub(c,a))
        if vectorops.dot(n,vectorops.sub(centroid,a)) > 0:
            b,c = c,b
            n = vectorops.mul(n,-1)
        try:
            n = vectorops.mul(n,1.0/vectorops.norm(n))
            glNormal3f(*n)
        except ZeroDivisionError:
            pass
        glBegin(GL_TRIANGLES)
        glTexCoord3f(uvwscale[0]*(a[0]-vmin[0]),uvwscale[1]*(a[1]-vmin[1]),uvwscale[2]*(a[2]-vmin[2]))
        glVertex3f(*a)
        glTexCoord3f(uvwscale[0]*(b[0]-vmin[0]),uvwscale[1]*(b[1]-vmin[1]),uvwscale[2]*(b[2]-vmin[2]))
        glVertex3f(*b)
        glTexCoord3f(uvwscale[0]*(c[0]-vmin[0]),uvwscale[1]*(c[1]-vmin[1]),uvwscale[2]*(c[2]-vmin[2]))
        glVertex3f(*c)
        glEnd()

class PrettyHullRenderer:
    """A faster scipy.spatial.ConvexHull renderer that uses caching and a nice rainbow effect."""
    def __init__(self,hull):
        self.hull = hull
        self.gl_object = glcommon.CachedGLObject()
        #initialized first render
        self.volumeTexture = None

    def render(self,transform=None,scale=1.0) :
        if self.volumeTexture is not None:
            glBindTexture(GL_TEXTURE_3D,self.volumeTexture[0])
        else:
            self.volumeTexture = glGenTextures(2)

            m,n,p = 1,512,1
            data = [0]*(m*n*p*4)

            glBindTexture(GL_TEXTURE_3D,self.volumeTexture[1])
            for i in range(n):
                data[i*4] = max(0,255-i)
                data[i*4+1] = 255-i if i >= 256 else i
                data[i*4+2] = max(i-255,0)
                data[i*4+3] = 255
            glTexImage3D(GL_TEXTURE_3D,0,GL_RGBA8,n,m,p,0,GL_RGBA,GL_UNSIGNED_BYTE,data)
            glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MAG_FILTER,GL_LINEAR)
            glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MIN_FILTER,GL_LINEAR)
            glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_S,GL_CLAMP);
            glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_T,GL_CLAMP);
            #glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_U,GL_CLAMP);

            glBindTexture(GL_TEXTURE_3D,self.volumeTexture[0])
            data = [0]*(m*n*p*4)
            for i in range(n):
                data[i*4+0] = 0
                data[i*4+1] = 0
                data[i*4+2] = 0
                data[i*4+3] = 0
                if (i+50) % 100 == 0:
                    data[i*4+0] = 0
                    data[i*4+1] = 0
                    data[i*4+2] = 0
                    data[i*4+3] = 255 * (512-i/2) / 512
            glTexImage3D(GL_TEXTURE_3D,0,GL_RGBA8,m,n,p,0,GL_RGBA,GL_UNSIGNED_BYTE,data)
            glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MAG_FILTER,GL_LINEAR)
            glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MIN_FILTER,GL_LINEAR)
            #glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_S,GL_REPEAT)
            #glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_T,GL_REPEAT)
            glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_S,GL_CLAMP);
            glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_T,GL_CLAMP);
            #glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_WRAP_U,GL_CLAMP);

        glEnable(GL_TEXTURE_3D)
        glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE)
        glDisable(GL_CULL_FACE)
        glDisable(GL_DEPTH_TEST)
        #glEnable(GL_LIGHTING)
        glDisable(GL_LIGHTING)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)
        #draw lines
        glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,(1,1,1,1))
        self.gl_object.draw(draw_hull,transform=transform,args=(self.hull,scale))
        glEnable(GL_CULL_FACE)
        glEnable(GL_DEPTH_TEST)
        #draw volume
        glEnable(GL_LIGHTING)
        glBindTexture(GL_TEXTURE_3D,self.volumeTexture[1])
        #glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,(1,0.5,0,0.8))
        glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,(1,1,1,0.8))
        self.gl_object.draw(draw_hull,transform=transform,args=(self.hull,scale))
        glDisable(GL_BLEND)
        glDisable(GL_TEXTURE_3D)

