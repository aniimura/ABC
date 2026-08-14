-- D9 fixture: Gap/ の structure に GapRecord(falsifier 込み)が無い(落ちるべき)
namespace Fixture
structure UnrecordedGap (n : Nat) : Prop where
  extra : n > 0
end Fixture
