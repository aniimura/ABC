import ABC3.Found.PGC.AdjoinFieldConstruction

/-!
# 代数閉包の同一視 `(K(x))‾ ≃ K‾`

`adjoinField K x` を `PAdicLocalField` として見ると、その代数閉包は
`AlgebraicClosure ↥K⟮x⟯` という**別の型**になる。しかし `K.closure` 自身が
`K(x)` 上の代数閉包でもあるので、両者は同型である。

これが要る理由: Lubin-Tate 塔は `K(Λ_n)/K(Λ_{n-1})` という**相対的な**拡大の
連なりとして与えられる(`ψ_n` は `𝒪_{K(Λ_{n-1})}` 上の Eisenstein)。
一方 `Found/PGC/TotallyRamified.lean` の完全分岐の道具は
`K.closure` の中の中間体について書かれている。両者を繋ぐには
`(K(Λ_{n-1}))‾ ≃ K‾` が要る。

## 配管の注意

`Algebra ↥K⟮x⟯ K.closure` などのインスタンスは、**中間体を直接書けば**
見つかるが `(adjoinField K x).carrier` と書くと見つからない
(`adjoinField` は `def` なので instance 探索が展開しない)。
そこで明示的に `instance` として与え直す。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped Valued

variable {p : ℕ} [Fact p.Prime]

/-- `K.closure` は `K(x)`-代数(中間体を直接書けば見つかるインスタンスの言い換え)。 -/
noncomputable instance algebraAdjoinFieldClosure (K : PAdicLocalField p) (x : K.closure) :
    Algebra ((adjoinField K x).carrier) K.closure :=
  (inferInstance : Algebra (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) K.closure)

instance isAlgebraicAdjoinFieldClosure (K : PAdicLocalField p) (x : K.closure) :
    Algebra.IsAlgebraic ((adjoinField K x).carrier) K.closure :=
  (inferInstance :
    Algebra.IsAlgebraic (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) K.closure)

/-- `K.closure` は `K(x)` 上の**代数閉包**でもある。 -/
instance isAlgClosureAdjoinField (K : PAdicLocalField p) (x : K.closure) :
    IsAlgClosure ((adjoinField K x).carrier) K.closure := by
  constructor
  · infer_instance
  · infer_instance

/-- **★★★★★代数閉包の同一視** `(K(x))‾ ≃ K‾`(`K(x)`-代数として)。

これで「`K(x)` 上の拡大」を `K.closure` の中の中間体として捉え直せる
——Lubin-Tate 塔の相対的な段を `K` からの絶対的な塔に繋ぐ鍵。 -/
noncomputable def closureEquivAdjoinField (K : PAdicLocalField p) (x : K.closure) :
    (adjoinField K x).closure ≃ₐ[(adjoinField K x).carrier] K.closure :=
  IsAlgClosure.equiv _ _ _

end ABC3.Found.PGC
