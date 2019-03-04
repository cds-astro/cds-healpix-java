package cds.healpix.common.sphgeom;


import static cds.healpix.Healpix.checkLatitude;
import static cds.healpix.common.math.Math.cos;
import static cds.healpix.common.math.Math.sin;
import static cds.healpix.common.math.Math.sqrt;
import static cds.healpix.common.math.Math.abs;
import static cds.healpix.common.math.Math.atan2;

/**
 * This class defines an elliptical cone.
 * The cone is defined by a center (ra, dec), semi-major and semi-minor axis (a, b) which are angular
 * distances (LOWER THAN pi/2!), and a position angle PA.
 * In the SIN (or  orthogonal) projection, centered around (ra=0, dec=0), i.e. (x=1, y=0, z=0), the
 * points (x, y, z) inside the ellise satifies:
 *  1 / (1 - rho^2) * [ x^2 / sig_x^2 - 2 * rho * x * y / (sig_x * sig_y)+  y^2 / sig_y^2 ]
 *  with:
 *  - sig_x^2 = sin^2(a)sin^2(PA) + sin^(b)cos^2(PA)
 *  - sig_y^2 = sin^2(a)cos^2(PA) + sin^(b)sin^2(PA)
 *  - p* sig_x * sigy = cos(PA)sin(PA)(sin^2(a) - sin^2(b)) 
 * 
 * 
 * 
 * FOR THE MOC: compute one value of circumcircle per order (like in the RUST lib)
 * => 2 ellipses per order (with +r_circul and -r_cirum)
 * => and use the RUST cone approx method, but for the elliptivcal cone :)
 * 
 * @author F.-X. Pineau
 *
 */
public class EllipticalCone {

  public static final double EPSILON = 1e-14;

  /** Longitude and lattitude of the center of elliptical cone. */
  protected double lon, lat;

  /** semi-major axis, semi-minor axis and position angle. */
  protected double a, b, pa;

  /** Components of the 3x3 rotation matrix transforming a vector into the
   * reference frame to a vector into the local frame (i.e. the frame in
   * which the position of the projection center is (1, 0, 0).
   * Remarks:
   *  r22 =  cos(lon)
   *  r21 = -sin(lon)
   *  r33 =  cos(lat)
   *  r13 =  sin(lat)
   */
  private double
  r11, r12, r13,
  r21, r22, r23,
  r31, r32, r33;

  private double sina, sinb;
  private double cosPA, sinPA;
  private double sigX2, sigY2, rho;
  private double oneOver1minRho2, twiceRhoOverSigXSigY;

  /** Component of the 2x2 rotation matrix transforming a point in the ellipse canonical frame
   * (in which rho = 0) to the ellipse frame rotated by the position angle. 
   * 
   */
  private double 
  a11, a12,
  a21, a22;
  

  public EllipticalCone(final double lonRad, final double latRad, 
      final double aRad, final double bRad) {
    this(lonRad, latRad, aRad, bRad, 0.0); 
  }

  public EllipticalCone(final double lonRad, final double latRad, 
      final double aRad, final double bRad, final double posAngRad) {
    this.setProjCenter(lonRad, latRad);
    this.setEllipseParams(aRad, bRad, posAngRad);
  }

  public double getA() {
    return this.a;
  }
  
  public double getB() {
    return this.b;
  }
  
  public double getSinA() {
    return this.sina;
  }
  
  public double getSinB() {
    return this.sinb;
  }
  
  /**
   * Returns {@code true} if the given point on the unit sphere is inside the elliptical cone.
   * @param lonRad longitude of the points, in radians
   * @param latRad latittude of the points, in radians
   * @return {@code true} if the given point on the unit sphere is inside the elliptical cone.
   */
  public boolean contains(final double lonRad, final double latRad) {
    final double[] projXY = new double[2]; 
    this.proj(lonRad, latRad, projXY);
    return squaredMahalanobisDistance(projXY[0], projXY[1]) <= 1;
  }
  
  /**
   * Returns {@code true} if the given cone overlap the elliptical cone.
   * @param lonRad longitude of the center of the cone, in radians
   * @param latRad latittude of the center of the cone, in radians
   * @param rRad cone radius, in radians
   * @return {@code true} if the given cone overlap the elliptical cone.
   */
  public boolean overlapCone(final double lonRad, final double latRad, final double rRad) {
    final double[] projXY = new double[2]; 
    this.proj(lonRad, latRad, projXY);
    return squaredMahalanobisDistance(projXY[0], projXY[1], rRad) <= 1;
  }
  
  /**
   * Returns {@code true} if the given cone is fully inside the elliptical cone.
   * @param lonRad longitude of the center of the cone, in radians
   * @param latRad latittude of the center of the cone, in radians
   * @param rRad cone radius, in radians
   * @return {@code true} if the given cone is fully inside the elliptical cone.
   */
  public boolean containsCone(final double lonRad, final double latRad, final double rRad) {
    final double[] projXY = new double[2]; 
    this.proj(lonRad, latRad, projXY);
    double d = squaredMahalanobisDistance(projXY[0], projXY[1], -rRad);
    return d == d && d <= 1;
  }
  
  /**
   * Returns the coordinates (lonRad, latRad) of points which are on the path along the ellipse edge
   * on the unit sphere. 
   * @param halfNumberOfPoints half the wanted number of points in the path
   * @return the coordinates (lonRad, latRad) of points which are on the path along the ellipse edge
   * on the unit sphere.
   */
  public double[][] pathAlongEdge(final int halfNumberOfPoints) {
    final double[][] coos = new double[2 * halfNumberOfPoints][];
    final double step = 2 * sina / halfNumberOfPoints;
    coos[0] = deproj(rotateEllipse(new double[]{sina, 0}));
    coos[halfNumberOfPoints] = deproj(rotateEllipse(new double[]{-sina, 0}));
    for (int i = 1; i < halfNumberOfPoints;) {
      double x = sina - i * step;
      x /= sina;
      double y = sinb * sqrt(1 - x * x);
      coos[i] = deproj(rotateEllipse(new double[]{x, y}));
      coos[coos.length - ++i] = deproj(rotateEllipse(new double[]{-x, -y}));
    }
    return coos;
  }
  
  private double[] rotateEllipse(final double[] xy) {
    double x = xy[0];
    double y = xy[1];
    xy[0] = this.a11 * x + this.a12 * y;
    xy[1] = this.a21 * x + this.a22 * y;
    return xy;
  }
  
  private double squaredMahalanobisDistance(final double x, final double y) {
    return oneOver1minRho2 * (
           (x * x) / this.sigX2 
        - twiceRhoOverSigXSigY * x * y
        +  (y * y) / this.sigY2
        );
  }
  
  private double squaredMahalanobisDistance(final double x, final double y, final double r) {
    if (r < 0 && (-r > this.a || -r > this.b)) {
      return Double.NaN;
    }
    double sinar = sin(this.a + r);
    double sinbr = sin(this.b + r);
    final double sa2 = sinar * sinar;
    final double sb2 = sinbr * sinbr;
    final double cpa2 = this.cosPA * this.cosPA;
    final double spa2 = this.sinPA * this.sinPA;
    double sigX2b = sa2 * spa2 + sb2 * cpa2;
    double sigY2b = sa2 * cpa2 + sb2 * spa2;
    double rhob = this.cosPA * this.sinPA * (sa2 - sb2);
    final double oneOver1minRho2b = 1 / (1 - rhob * rhob);
    final double twiceRhoOverSigXSigYb = 2 * rhob / sqrt(sigX2b * sigY2b);
    return oneOver1minRho2b * (
        (x * x) / sigX2b
     - twiceRhoOverSigXSigYb * x * y
     +  (y * y) / sigY2b
    );
  }
  
  public void setProjCenter(final double lon, final double lat) {
    checkLatitude(lat);
    this.lon = lon;
    this.lat = lat;
    final double ca = cos(lon); // ca stands for cos(alpha)
    final double sa = sin(lon); // sa stands for sin(alpha)
    final double cd = cos(lat); // cd stands for cos(delta)
    final double sd = sin(lat); // sd stands for sin(delta)
    this.r11 =  ca * cd; this.r12 =  sa * cd; this.r13 = sd;
    this.r21 =      -sa; this.r22 =       ca; this.r23 =  0;
    this.r31 = -ca * sd; this.r32 = -sa * sd; this.r33 = cd;
  }


  private void setEllipseParams(final double a, final double b, final double pa) {
    assert a > 0 && b > 0;
    this.a = a;
    this.b = b;
    this.pa = pa;
    this.sina = sin(a);
    final double sa2 = this.sina * this.sina;
    this.sinb = sin(b);
    final double sb2 = this.sinb * this.sinb;
    this.cosPA = cos(pa);
    this.sinPA =  sin(pa);
    final double cpa2 = this.cosPA * this.cosPA;
    final double spa2 = this.sinPA * this.sinPA;
    this.sigX2 = sa2 * spa2 + sb2 * cpa2;
    this.sigY2 = sa2 * cpa2 + sb2 * spa2;
    this.rho = this.cosPA * this.sinPA * (sa2 - sb2);
    this.oneOver1minRho2 = 1 / (1 - this.rho * this.rho);
    this.twiceRhoOverSigXSigY = 2 * this.rho / sqrt(this.sigX2 * this.sigY2);
    
    this.a11 = this.sinPA; this.a12 = -this.cosPA;    
    this.a21 = this.cosPA; this.a22 =  this.sinPA;
  }

  private boolean proj(final double lon, final double lat, final double[] resultXY) {
    checkLatitude(lat);
    // Transforms LonLat into a 3d vector
    double z = cos(lat);
    double x = z * cos(lon);
    double y = z * sin(lon);
    z = sin(lat);
    // Assertions to verify coordinate are coherent with expected values
    assert -1 <= x && x <= 1 : "x: " + x + " should be in [-1, 1]";
    assert -1 <= y && y <= 1 : "y: " + y + " should be in [-1, 1]";
    assert -1 <= z && z <= 1 : "z: " + z + " should be in [-1, 1]";
    assert abs(1 - (x * x + y * y + z * z)) < EPSILON : "(x^2 + y^2 + z^2): " + (x * x + y * y + z * z);
    // Rotate the vector to the local frame and performs the projection
    return this.projAssert(
        this.r11 * x + this.r12 * y + this.r13 * z,
        this.r21 * x + this.r22 * y + this.r23 * z,
        this.r31 * x + this.r32 * y + this.r33 * z,
        resultXY);
  }

  private final boolean projAssert(final double x, final double y, final double z, final double[] resultXY) {
    assert -1 <= x && x <= 1 : "x: " + x + " should be in [-1, 1]";
    assert -1 <= y && y <= 1 : "y: " + y + " should be in [-1, 1]";
    assert -1 <= z && z <= 1 : "z: " + z + " should be in [-1, 1]";
    assert abs(1 - (x * x + y * y + z * z)) < EPSILON : "(x^2 + y^2 + z^2): " + (x * x + y * y + z * z);
    return this.proj(x, y, z, resultXY);
  }

  protected boolean proj(final double x, final double y, final double z, final double[] resultXY) {
    if (x < 0) { // Back hemisphere
      resultXY[0] = Double.NaN;
      resultXY[1] = Double.NaN;
      return false;
    }
    resultXY[0] = y;
    resultXY[1] = z;
    return true;
  }

  private double[] deproj(final double[] xy2lonlat) {
    final double x2D = xy2lonlat[0];
    final double y2D = xy2lonlat[1];
    final double r2 = x2D * x2D + y2D * y2D;
    double x2 = 1 - r2; // = x^2 + y^2 + z^2 = 1
    if (x2 < 0) { // Back hemisphere 
      if (x2 > -2e-16) { // But accept some rounding error
        x2 = 0;
      } else {
        xy2lonlat[0] = Double.NaN;
        xy2lonlat[0] = Double.NaN;
        return xy2lonlat;
      }
    }
    this.setLonLatFromLocalXYZ(xy2lonlat, sqrt(x2), x2D, y2D);
    return xy2lonlat;
  }

  /**
   * To be called at the end of the {@link #deproj(XY, SetableLonLat)}
   * method when local coordinates (x, y, z) have been computed.
   * @param toLonLat
   * @param x x coordinate in the local frame
   * @param y y coordinate in tlongitudeRadianshe local frame
   * @param z z coordinate in the local frame
   */
  private final void setLonLatFromLocalXYZ(final double[] toLonLat,
      final double x, final double y, final double z) {
    // Rotate the vector to the global frame (using the transpose of R)
    setLonLatFromXYZ(toLonLat,
        this.r11 * x + this.r21 * y + this.r31 * z,
        this.r12 * x + this.r22 * y + this.r32 * z,
        this.r13 * x + this.r23 * y + this.r33 * z);
  }

  private static final void setLonLatFromXYZ(final double[] toLonLat,
      final double x, final double y, final double z) {
    assert -1 - EPSILON <= x && x <= 1 + EPSILON : "x: " + x + " should be in [-1, 1]";
    assert -1 - EPSILON <= y && y <= 1 + EPSILON : "y: " + y + " should be in [-1, 1]";
    assert -1 - EPSILON <= z && z <= 1 + EPSILON : "z: " + z + " should be in [-1, 1]";
    assert abs(1 - sqrt(x * x + y * y + z * z)) < EPSILON : "sqrt(x^2 + y^2 + z^2): " + sqrt(x * x + y * y + z * z);
    // Length of the projection on the xy plane
    double r2 = x * x + y * y;
    // Latitude in [-pi/2, pi/2] (ok, since cos always positive here)
    toLonLat[1] = atan2(z, sqrt(r2));
    // Compute the longitude in [-pi, pi]
    r2 = atan2(y, x);
    // Conforms to convention: Longitude in [0, 2*PI]
    toLonLat[0] = r2 < 0 ? 2 * Math.PI + r2 : r2;
  }

}
