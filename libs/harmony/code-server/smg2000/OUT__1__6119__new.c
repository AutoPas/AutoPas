void OUT 1 6755 ( void  o u t a r g v )
{
  int Ai = ( ( int ) ( o u t a r g v [ 2 0 ] ) ) ;
  int xi = ( ( int ) ( o u t a r g v [ 1 9 ] ) ) ;
  int r i = ( ( int ) ( o u t a r g v [ 1 8 ] ) ) ;
  double Ap = ( ( double ) ( o u t a r g v [ 1 7 ] ) ) ;
  double xp = ( ( double ) ( o u t a r g v [ 1 6 ] ) ) ;
  double rp = ( ( double ) ( o u t a r g v [ 1 5 ] ) ) ;
  int l o o p i = ( ( int ) ( o u t a r g v [ 1 4 ] ) ) ;
  int l o o p j = ( ( int ) ( o u t a r g v [ 1 3 ] ) ) ;
  int loopk = ( ( int ) ( o u t a r g v [ 1 2 ] ) ) ;
  int hypr e sx1 = ( ( int ) ( o u t a r g v [ 1 1 ] ) ) ;
  int hypr e sy1 = ( ( int ) ( o u t a r g v [ 1 0 ] ) ) ;
  int hypr e s z 1 = ( ( int ) ( o u t a r g v [ 9 ] ) ) ;
  int hypr e sx2 = ( ( int ) ( o u t a r g v [ 8 ] ) ) ;
  int hypr e sy2 = ( ( int ) ( o u t a r g v [ 7 ] ) ) ;
  int hypr e s z 2 = ( ( int ) ( o u t a r g v [ 6 ] ) ) ;
  int hypr e sx3 = ( ( int ) ( o u t a r g v [ 5 ] ) ) ;
  int hypr e sy3 = ( ( int ) ( o u t a r g v [ 4 ] ) ) ;
  int hypr e s z 3 = ( ( int ) ( o u t a r g v [ 3 ] ) ) ;
  int hypr e nx = ( ( int ) ( o u t a r g v [ 2 ] ) ) ;
  int hypr e ny = ( ( int ) ( o u t a r g v [ 1 ] ) ) ;
  int hypr e nz = ( ( int ) ( o u t a r g v [ 0 ] ) ) ;
  for ( loopk = 0 ; loopk < hypr e nz ; loopk++) {
    for ( l o o p j = 0 ; l o o p j < hypr e ny ; l o o p j++) {
      for ( l o o p i = 0 ; l o o p i < hypr e nx ; l o o p i++) {{
	rp [ r i ] â= ( (Ap[ Ai ] )  ( xp [ xi ] ) ) ;
      }
      Ai += hypr e sx1 ;
      xi += hypr e sx2 ;
      r i += hypr e sx3 ;
      }
      Ai += ( hypr e sy1 â ( hypr e nx  hypr e sx1 ) ) ;
      xi += ( hypr e sy2 â ( hypr e nx  hypr e sx2 ) ) ;
      r i += ( hypr e sy3 â ( hypr e nx  hypr e sx3 ) ) ;
    }
    Ai += ( hypr e s z 1 â ( hypr e ny  hypr e sy1 ) ) ;
    xi += ( hypr e s z 2 â ( hypr e ny  hypr e sy2 ) ) ;
    r i += ( hypr e s z 3 â ( hypr e ny  hypr e sy3 ) ) ;
  }
  ( ( int ) ( o u t a r g v [ 0 ] ) ) = hypr e nz ;
  ( ( int ) ( o u t a r g v [ 1 ] ) ) = hypr e ny ;
  ( ( int ) ( o u t a r g v [ 2 ] ) ) = hypr e nx ;
  ( ( int ) ( o u t a r g v [ 3 ] ) ) = hypr e s z 3 ;
  ( ( int ) ( o u t a r g v [ 4 ] ) ) = hypr e sy3 ;
  ( ( int ) ( o u t a r g v [ 5 ] ) ) = hypr e sx3 ;
  ( ( int ) ( o u t a r g v [ 6 ] ) ) = hypr e s z 2 ;
  ( ( int ) ( o u t a r g v [ 7 ] ) ) = hypr e sy2 ;
  ( ( int ) ( o u t a r g v [ 8 ] ) ) = hypr e sx2 ;
  ( ( int ) ( o u t a r g v [ 9 ] ) ) = hypr e s z 1 ;
  ( ( int ) ( o u t a r g v [ 1 0 ] ) ) = hypr e sy1 ;
  ( ( int ) ( o u t a r g v [ 1 1 ] ) ) = hypr e sx1 ;
  ( ( int ) ( o u t a r g v [ 1 2 ] ) ) = loopk ;
  ( ( int ) ( o u t a r g v [ 1 3 ] ) ) = l o o p j ;
  ( ( int ) ( o u t a r g v [ 1 4 ] ) ) = l o o p i ;
  ( ( double ) ( o u t a r g v [ 1 5 ] ) ) = rp ;
  ( ( double ) ( o u t a r g v [ 1 6 ] ) ) = xp ;
  ( ( double ) ( o u t a r g v [ 1 7 ] ) ) = Ap;
  ( ( int ) ( o u t a r g v [ 1 8 ] ) ) = r i ;
  ( ( int ) ( o u t a r g v [ 1 9 ] ) ) = xi ;
  ( ( int ) ( o u t a r g v [ 2 0 ] ) ) = Ai ;
}	
