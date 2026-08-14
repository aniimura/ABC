-- D32 fixture: 文字列を潰しても**本物の** sorry は落とせる(落ちるべき)
-- ★D31 と対にする。「全部通す」器具にならないことの確認。
namespace Fixture
theorem realSorry : True := by
  have h : True := sorry
  exact h
end Fixture
