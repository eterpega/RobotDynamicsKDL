classdef SpatialForceVectorSemanticsBase < SwigRef
  methods
    function this = swig_this(self)
      this = iDynTreeMEX(3, self);
    end
    function self = SpatialForceVectorSemanticsBase(varargin)
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = iDynTreeMEX(524, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.SwigClear();
      end
    end
    function varargout = check_linear2angularConsistency(self,varargin)
      [varargout{1:nargout}] = iDynTreeMEX(525, self, varargin{:});
    end
    function varargout = toString(self,varargin)
      [varargout{1:nargout}] = iDynTreeMEX(526, self, varargin{:});
    end
    function varargout = display(self,varargin)
      [varargout{1:nargout}] = iDynTreeMEX(527, self, varargin{:});
    end
    function delete(self)
      if self.swigPtr
        iDynTreeMEX(528, self);
        self.SwigClear();
      end
    end
  end
  methods(Static)
  end
end
