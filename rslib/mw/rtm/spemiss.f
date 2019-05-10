      subroutine spemiss(freq,tk,theta,sal,ev,eh)
      
      real freq,theta,tk,sal,wind,ev,eh
      
      wind = 0.0
      
      call rough_ks(freq,theta,tk,sal,wind,ev,eh)
      
      end
