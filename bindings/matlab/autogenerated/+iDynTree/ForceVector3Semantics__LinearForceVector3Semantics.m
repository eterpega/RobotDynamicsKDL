classdef ForceVector3Semantics__LinearForceVector3Semantics < iDynTree.GeomVector3Semantics__LinearForceVector3Semantics
  methods
    function self = ForceVector3Semantics__LinearForceVector3Semantics(varargin)
      self@iDynTree.GeomVector3Semantics__LinearForceVector3Semantics(SwigRef.Null);
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = iDynTreeMEX(478, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.SwigClear();
      end
    end
    function delete(self)
      if self.swigPtr
        iDynTreeMEX(481, self);
        self.SwigClear();
      end
    end
  end
  methods(Static)
    function varargout = compose(varargin)
     [varargout{1:nargout}] = iDynTreeMEX(479, varargin{:});
    end
    function varargout = inverse(varargin)
     [varargout{1:nargout}] = iDynTreeMEX(480, varargin{:});
    end
  end
end
