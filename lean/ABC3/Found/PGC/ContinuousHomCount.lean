import Mathlib.Topology.Algebra.ContinuousMonoidHom
import Mathlib.Topology.Algebra.OpenSubgroup
import Mathlib.Data.ZMod.Basic

/-!
# 連続準同型の個数 `N_n(G) = #Hom_cont(G, ℤ/n)`

[pGC] Proposition 1.2 への経路 C(`ResearchPaper/pgc-goal.md`)は、
絶対 Galois 群 `Γ_K` の**位相群としての同型類だけ**から `q` と `[K:ℚ_p]` を
読み取ることを目指す。そのための不変量が

  `N_n(G) := # Hom_cont(G, ℤ/n)`

である。本ファイルはこの不変量を定義し、位相群の同型で不変であることを示す。

## `ContinuousMonoidHom` を使わない理由(逸脱の記録)

mathlib の `ContinuousMonoidHom G A` は `A` に位相を要求するが、
`ZMod n` には `TopologicalSpace` インスタンスが**無い**。離散位相を
`letI` で持ち込むと、`Multiplicative (ZMod n)` 上の同型を経由するたびに
インスタンスの一致を証明する羽目になる。

そこで「連続」を**核が開である**ことで定義する:

  `f : G →* A` が連続  ⟺  `IsOpen (ker f)`.

`A` が離散位相を持つときこれは通常の連続性と同値であり、しかも `A` の
位相を一切参照しない。指標(`A = ℤ/n`)を扱う本経路ではこれで十分である。

## `SeparatelyContinuousMul G` を仮定する理由(逸脱の記録)

`ker (f * g) ⊇ ker f ⊓ ker g` から `ker (f * g)` の開性を出すには
「開部分群を含む部分群は開」(`Subgroup.isOpen_mono`)が要る。mathlib の
その補題は `[SeparatelyContinuousMul G]` を仮定しているので、こちらも
同じ仮定を置く。`IsTopologicalGroup G` から自動で得られるので、
`Γ_K`(Krull 位相)への適用には何の障害も無い。
-/

namespace ABC3.Found.PGC

/-- 連続(= 核が開)な準同型のなす `G →* A` の部分群。 -/
def contHom (G : Type*) [Group G] [TopologicalSpace G] [SeparatelyContinuousMul G]
    (A : Type*) [CommGroup A] : Subgroup (G →* A) where
  carrier := {f | IsOpen ((MonoidHom.ker f : Subgroup G) : Set G)}
  one_mem' := by
    simp only [Set.mem_setOf_eq, MonoidHom.ker_one]
    simp
  mul_mem' := by
    intro f g hf hg
    simp only [Set.mem_setOf_eq] at hf hg ⊢
    refine Subgroup.isOpen_mono (H₁ := MonoidHom.ker f ⊓ MonoidHom.ker g) ?_ ?_
    · intro x hx
      simp only [Subgroup.mem_inf, MonoidHom.mem_ker] at hx ⊢
      simp [hx.1, hx.2]
    · simpa using hf.inter hg
  inv_mem' := by
    intro f hf
    simp only [Set.mem_setOf_eq] at hf ⊢
    convert hf using 2
    ext x
    simp [MonoidHom.mem_ker]

@[simp] lemma mem_contHom {G : Type*} [Group G] [TopologicalSpace G] [SeparatelyContinuousMul G]
    {A : Type*} [CommGroup A] (f : G →* A) :
    f ∈ contHom G A ↔ IsOpen ((MonoidHom.ker f : Subgroup G) : Set G) := Iff.rfl

/-- 連続準同型 `φ` を右から合成しても連続性は保たれる(核の引き戻しが開)。 -/
lemma comp_mem_contHom {G H : Type*} [Group G] [TopologicalSpace G] [SeparatelyContinuousMul G]
    [Group H] [TopologicalSpace H] [SeparatelyContinuousMul H] {A : Type*} [CommGroup A]
    (φ : G →* H) (hφ : Continuous φ) {f : H →* A} (hf : f ∈ contHom H A) :
    f.comp φ ∈ contHom G A := by
  rw [mem_contHom] at hf ⊢
  rw [← MonoidHom.comap_ker, Subgroup.coe_comap]
  exact hf.preimage hφ

/-- 位相群の同型 `α : G ≃ₜ* H` に沿った連続準同型の対応 `f ↦ f ∘ α`。 -/
def contHomEquiv {G H : Type*} [Group G] [TopologicalSpace G] [SeparatelyContinuousMul G]
    [Group H] [TopologicalSpace H] [SeparatelyContinuousMul H] (α : G ≃ₜ* H)
    (A : Type*) [CommGroup A] : contHom H A ≃ contHom G A where
  toFun f := ⟨(f : H →* A).comp α.toMulEquiv.toMonoidHom,
    comp_mem_contHom _ (map_continuous α) f.2⟩
  invFun g := ⟨(g : G →* A).comp α.symm.toMulEquiv.toMonoidHom,
    comp_mem_contHom _ (map_continuous α.symm) g.2⟩
  left_inv f := by
    refine Subtype.ext (MonoidHom.ext fun x => ?_)
    show (f : H →* A) (α.toMulEquiv (α.toMulEquiv.symm x)) = (f : H →* A) x
    simp
  right_inv g := by
    refine Subtype.ext (MonoidHom.ext fun x => ?_)
    show (g : G →* A) (α.toMulEquiv.symm (α.toMulEquiv x)) = (g : G →* A) x
    simp

/-- 不変量 `N_n(G) = # Hom_cont(G, ℤ/n)`。 -/
noncomputable def contHomCard (G : Type*) [Group G] [TopologicalSpace G]
    [SeparatelyContinuousMul G] (n : ℕ) : ℕ :=
  Nat.card (contHom G (Multiplicative (ZMod n)))

/-- `N_n` は位相群の同型で不変。 -/
theorem contHomCard_congr {G H : Type*} [Group G] [TopologicalSpace G] [SeparatelyContinuousMul G]
    [Group H] [TopologicalSpace H] [SeparatelyContinuousMul H] (α : G ≃ₜ* H) (n : ℕ) :
    contHomCard G n = contHomCard H n :=
  (Nat.card_congr (contHomEquiv α (Multiplicative (ZMod n)))).symm

end ABC3.Found.PGC
