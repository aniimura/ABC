/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.DivisorOfSectionEq
import ABC3.Found.GenEll.GlobalRatio
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★段 E0 の要 —— チャートの上で `div(s₀)` は `(s₀/t)` である（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★★★★★★★★★★★これは何か —— 段 E0 の壁が消えた

段 E0（切断の零点因子 `divisorOfSection`）の障害は
「`divIdeal` が**自明化の無い**アフィン開で `⊤` になること」であった（`§9-809`）。
★`§9-884` は「大域的に自明なら等号」まで来ていたが、
一般の局所自明な `M` については**`X_{s_i}` 上の自明化**が要ると見えていた。

★★★★**しかしそれは要らなかった**（2026-08-29 実測）:

    `V` が自明化つきアフィン開なら、`X_t ⊓ V` の上で
    **`divIdeal M s₀ = span { (s₀/t) の制限 }`**

——`trivValue(s₀)` と `sectionRatio(s₀, t)` は**単元 `trivValue(t)⁻¹` 倍しか違わない**ので、
生成するイデアルは同じである。★`t` で正規化した自明化を**作る必要が無い**。

## ★★★これで何が繋がるか

★`X_t ⊓ V = X.basicOpen (trivValue M V e t)` （`nonVanishing_inf`、段 D2）なので、
**`V` がアフィンなら `X_t ⊓ V` もアフィン**（基本開集合）である。
★★したがって `{X_{s_i} ⊓ U_j}` は

* `X` の**アフィン開被覆**（`⨆_i X_{s_i} = ⊤` かつ `⨆_j U_j = ⊤`）
* どの成分にも**自明化がある**
* どの成分でも `divIdeal M s₀ = span {s₀/s_i}`

を同時に満たす。★★★mathlib の `Scheme.IdealSheafData.ext_of_iSup_eq_top`
（アフィン被覆の上で一致すれば等しい）に渡せる形である。

## ★残っている段（明示）

★★★★残るのは **`(ψ^*超平面).ideal` を同じ被覆の上で計算すること**である。
`§9-921` は**点ごと**（`pullbackIdeal`）には `span {(s₀/s_i)(x)}` だと言っているので、
イデアル層の水準へ上げれば `divisorOfSection M s₀ = ψ^*超平面` が出る。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★★★自明化つきアフィン開との交わりはアフィン -/

/-- ★★★**`X_t ⊓ V` はアフィンである**（`V` がアフィンで自明化を持てば）。

★`nonVanishing_inf`（段 D2）で `X.basicOpen (trivValue M V e t)` に書き換わるので、
アフィン開の基本開集合としてアフィンである。 -/
theorem isAffineOpen_nonVanishing_inf (M : X.PresheafOfModules)
    (V : X.Opens) (hV : IsAffineOpen V)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (t : (M.obj (op ⊤) : Type)) :
    IsAffineOpen (nonVanishing M t ⊓ V) := by
  rw [nonVanishing_inf M V e t]
  exact hV.basicOpen _

/-- ★★**`{X_{s_i} ⊓ U_j}` は `X` を覆う**。 -/
theorem iSup_nonVanishing_inf_eq_top {ι κ : Type} (M : X.PresheafOfModules)
    (s : κ → (M.obj (op ⊤) : Type)) (U : ι → X.Opens)
    (hcovU : (⨆ j, U j) = ⊤) (hcovs : (⨆ i, nonVanishing M (s i)) = ⊤) :
    (⨆ p : κ × ι, nonVanishing M (s p.1) ⊓ U p.2) = ⊤ := by
  refine eq_top_iff.2 ?_
  intro x _
  have h1 : x ∈ (⨆ i, nonVanishing M (s i)) := by rw [hcovs]; trivial
  have h2 : x ∈ (⨆ j, U j) := by rw [hcovU]; trivial
  obtain ⟨i, hi⟩ := TopologicalSpace.Opens.mem_iSup.1 h1
  obtain ⟨j, hj⟩ := TopologicalSpace.Opens.mem_iSup.1 h2
  exact TopologicalSpace.Opens.mem_iSup.2 ⟨(i, j), ⟨hi, hj⟩⟩

/-! ## ★★★★★★★★★★★★★段 E0 の要 -/

/-- ★★★★★★★★★★★★★**`X_t ⊓ V` の上で `divIdeal M s₀` は `(s₀/t)` が生成する**。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★★★★**要点**——`t` で正規化した自明化を**作る必要が無い**。
`trivValue(s₀)` と `sectionRatio(s₀, t)` は**単元 `trivValue(t)⁻¹` 倍しか違わない**ので、
生成するイデアルは同じだからである。

★これで段 E0 の「自明化の無い開で `⊤` になる」障害が、
**自明化つきアフィン開の被覆の上では消える**。 -/
theorem divIdeal_nonVanishing_inf (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (V : X.Opens) (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s₀ t : (M.obj (op ⊤) : Type))
    (haff : IsAffineOpen (nonVanishing M t ⊓ V)) :
    divIdeal M s₀ ⟨nonVanishing M t ⊓ V, haff⟩
      = Ideal.span {X.presheaf.map (homOfLE (inf_le_left :
          nonVanishing M t ⊓ V ≤ nonVanishing M t)).op (globalRatio M hM s₀ t)} := by
  rw [divIdeal_eq M s₀ ⟨_, haff⟩
      (trivialOfLe M (inf_le_right : nonVanishing M t ⊓ V ≤ V) e),
    globalRatio_res M hM s₀ t ⟨V, e⟩,
    trivValue_restrict M (inf_le_right : nonVanishing M t ⊓ V ≤ V) e s₀,
    sectionRatio]
  exact (Ideal.span_singleton_mul_right_unit
    (isUnit_trivValue_res M V e t).unit⁻¹.isUnit _).symm

/-! ## ★出典の紐付け(`.src`) -/

def isAffineOpen_nonVanishing_inf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2(X_t ⊓ V はアフィンである)",
    sectionId := "genell-def-1-2" }

def iSup_nonVanishing_inf_eq_top.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2({X_{s_i} ⊓ U_j} は X を覆う)",
    sectionId := "genell-def-1-2" }

def divIdeal_nonVanishing_inf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2(X_t ⊓ V の上で divIdeal は (s₀/t) が生成する)",
    sectionId := "genell-def-1-2" }

def divIdeal_nonVanishing_inf.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "divIdeal_eq(自明化があればどれで測ってもよい、§9-809)"
      (.inProject "ABC3" "ABC3.Found.GenEll.divIdeal_eq") 1,
    .citation "[ABC3]" "globalRatio_res(大域の比はチャートで sectionRatio に戻る、§9-841)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalRatio_res") 1,
    .citation "[ABC3]" "trivValue_restrict(自明化の制限、段 D2)"
      (.inProject "ABC3" "ABC3.Found.GenEll.trivValue_restrict") 1,
    .implicitStep
      ("★★★★測定(2026-08-29): 段 E0 の障害は「divIdeal が自明化の無いアフィン開で ⊤ になる」" ++
       "ことであり、一般の局所自明な M については X_{s_i} 上の自明化を" ++
       "**作る**必要があると見えていた(それは前層の同型の貼り合わせで、枠組みに無い)。" ++
       "★**しかしそれは要らなかった**——trivValue(s₀) と sectionRatio(s₀, t) は" ++
       "**単元 trivValue(t)⁻¹ 倍しか違わない**ので、生成するイデアルは同じである") 5,
    .implicitStep
      ("★★これで {X_{s_i} ⊓ U_j} が『アフィン開被覆・自明化つき・divIdeal が既知』を" ++
       "同時に満たし、mathlib の Scheme.IdealSheafData.ext_of_iSup_eq_top に渡せる。" ++
       "★残るのは (ψ^*超平面).ideal を同じ被覆の上で計算することである") 4 ]

end ABC3.Found.GenEll
