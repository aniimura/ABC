import ABC3.Found.GenEll.ProjTopology

/-!
# [GenEll] Definition 1.1 の足場 —— **射影多様体はコンパクト**(GAGA なし)(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> There is an evident notion of tensor product of arithmetic line bundles on X. The

## ★★★`X^arc` のコンパクト性に GAGA は要らない

`Proposition 1.6` のアルキメデス側は
「コンパクト空間 `X^arc` 上の連続関数は有界」を使う(`ArchBound.lean`)。
★その `X^arc` のコンパクト性が要る。

`ProjTopology.lean` で `ℙⁿ(ℂ)` のコンパクト性は取った。
残るのは **`X(ℂ)` が `ℙⁿ(ℂ)` の閉集合である**ことである。

★★★**これは初等的に出る。** `X` は斉次多項式の零点集合であり、
**多項式は連続**だからその零点集合は閉である。
★★複素解析空間も GAGA も要らない——**連続性だけ**である。

## ★★機構は商位相そのもの

`ℙ(V)` は `{v ≠ 0}` の商である。商位相では
**閉集合 ⟺ 引き戻しが閉**(`isClosed_coinduced`)。

`V` の閉錐(スカラー倍で閉じた閉集合)`S` に対し、
`{p : ℙ(V) | p.rep ∈ S}` の引き戻しは `{v | ↑v ∈ S}` であり、
★`S` が錐なので代表元の取り方に依らない(`exists_smul_eq_mk_rep`)。

★★★**斉次多項式の零点集合は閉錐である**——これで話が閉じる。

## ★これで「正面から要ると思ったものが要らなかった」がまた 1 つ

複素解析空間 → 位相とコンパクト性だけ(`ProjTopology.lean`)に続いて、
**GAGA → 多項式の連続性だけ**である。
-/

namespace ABC3.Found.GenEll

open scoped LinearAlgebra.Projectivization

variable {V : Type*} [NormedAddCommGroup V] [NormedSpace ℂ V]

/-! ## ★★閉錐の像は閉 -/

/-- ★**錐**(原点を除いてスカラー倍で閉じた集合)。 -/
def IsCone (S : Set V) : Prop :=
  ∀ (c : ℂˣ) {v : V}, v ∈ S → (c : ℂ) • v ∈ S

/-- ★★★**閉錐が `ℙ(V)` に定める集合は閉**。

★機構は商位相そのもの(`isClosed_coinduced`)。
★★錐であることが「代表元の取り方に依らない」を保証する。 -/
theorem isClosed_projSetOfCone {S : Set V} (hS : IsClosed S) (hcone : IsCone S) :
    IsClosed {p : ℙ ℂ V | p.rep ∈ S} := by
  have hq : Topology.IsQuotientMap
      (@Quotient.mk' {w : V // w ≠ 0} (projectivizationSetoid ℂ V)) :=
    isQuotientMap_quotient_mk'
  refine hq.isClosed_preimage.mp ?_
  -- 引き戻しは `{v | ↑v ∈ S}`(錐なので代表元の取り方に依らない)
  convert hS.preimage (continuous_subtype_val (p := fun w : V => w ≠ 0)) using 1
  ext v
  obtain ⟨a, ha⟩ := Projectivization.exists_smul_eq_mk_rep ℂ (v : V) v.2
  constructor
  · intro hv
    -- `a • v = rep ∈ S` から `v ∈ S`(`a⁻¹` を掛ける)
    have h1 : (a : ℂ) • (v : V) ∈ S := ha ▸ hv
    have h2 : ((a⁻¹ : ℂˣ) : ℂ) • ((a : ℂ) • (v : V)) ∈ S := hcone a⁻¹ h1
    rwa [smul_smul, ← Units.val_mul, inv_mul_cancel, Units.val_one, one_smul] at h2
  · intro hv
    show (Projectivization.mk ℂ (v : V) v.2).rep ∈ S
    exact ha ▸ hcone a hv

/-! ## ★★★コンパクト性 -/

variable [FiniteDimensional ℂ V]

/-- ★★★**閉錐が定める射影集合はコンパクト**。

原文 (GenEll p.3):
> There is an evident notion of tensor product of arithmetic line bundles on X. The

★★`ℙ(V)` がコンパクト(`ProjTopology.lean`)で、その閉部分集合だからである。
★★★**これが `X^arc` のコンパクト性である**——複素解析空間も GAGA も要らない。 -/
theorem isCompact_projSetOfCone {S : Set V} (hS : IsClosed S) (hcone : IsCone S) :
    IsCompact {p : ℙ ℂ V | p.rep ∈ S} := by
  haveI := compactSpace_projectivization (V := V)
  exact (isClosed_projSetOfCone hS hcone).isCompact

omit [FiniteDimensional ℂ V] in
/-- ★★**零点集合は閉錐である**(連続かつ斉次な関数について)。

★これが「斉次多項式の零点集合」の抽象化である——
**多項式であることは使わず、連続性と斉次性だけを使う**。 -/
theorem isCone_zeroSet {f : V → ℂ} (hhom : ∀ (c : ℂˣ) (v : V), f ((c : ℂ) • v) = 0 ↔ f v = 0) :
    IsCone {v : V | f v = 0} := by
  intro c v hv
  exact (hhom c v).2 hv

/-- ★★★**連続かつ斉次な関数の零点集合はコンパクトな射影集合を定める**。

★★★**これが `X^arc` のコンパクト性の最終形**である。
`X ⊆ ℙⁿ` を切り出す斉次多項式に当てればよい——多項式は連続である。 -/
theorem isCompact_projZeroSet {f : V → ℂ} (hf : Continuous f)
    (hhom : ∀ (c : ℂˣ) (v : V), f ((c : ℂ) • v) = 0 ↔ f v = 0) :
    IsCompact {p : ℙ ℂ V | f p.rep = 0} :=
  isCompact_projSetOfCone (isClosed_eq hf continuous_const) (isCone_zeroSet hhom)

/-! ## ★出典の紐付け(`.src`) -/

def isCompact_projZeroSet.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1(X^arc のコンパクト性——射影集合の場合)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.GenEll
