Pod::Spec.new do |s|
  s.name         = "libapngiter"
  s.version      = "0.1"
  s.summary      = "iteration port libapng from maxst.users.sourceforge.net"
  s.homepage     = "https://github.com/maximgavrilov/libapngiter"
  s.license      = 'zlib'
  s.author       = { "Maxim Gavrilov" => "maxim.gavrilov@gmail.com" }
  s.source       = { :git => "https://github.com/maximgavrilov/libapngiter.git" }
  s.platform     = :ios
  s.source_files = 'libapngiter'
end
