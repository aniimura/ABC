import ABC3.Meta.Claim
import Mathlib.Algebra.Category.Ring.Basic
import Mathlib.CategoryTheory.Limits.Shapes.Pullback.IsPullback.Basic
import Mathlib.RingTheory.Ideal.Quotient.Operations

/-!
# Arakelov (B2) 第 194 ブロック —— **押し出しの核はイデアルの拡大**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★残り 2 欄(`comap`)の核心

§9-211 で測った通り、`isCartierDivisor_comap` と `ofDivisor_pullback` に要るのは

    (D.comap f).ideal A = (D.ideal B).map (f.appLE B A e)

であり、`Hom.ker_apply` で左辺は `RingHom.ker ((pullback.fst f D.subschemeι).app A)` に落ちる。
★★あとは**アフィンの引き戻しが環の押し出し**であること
(mathlib `isPushout_appTop_of_isPullback`)を使い、
本ブロックの純環論の補題で核を計算すればよい。

## ★★普遍性だけで出る

    A --u--> B
    |        |
    v        iB
    ↓        ↓
    C --iC--> P     （押し出し、v は全射で核 I）

    ⟹  ker iB = I.map u

| 向き | 論法 |
|---|---|
| ⊇ | `iB(u a) = iC(v a) = iC 0 = 0` ✅ 図式の可換性だけ |
| ⊆ | ★`B/(I·B)` への余錐を作り、普遍性で `P → B/(I·B)` を得る |

★★★テンソル積を一度も使わずに済む——`C ≅ A/I`(全射だから)と
`Ideal.Quotient.lift` で余錐が組めるからである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `ker_of_isPushout` | ★★★★**押し出しで片脚が全射なら反対の脚の核は拡大** |
-/

universe u

namespace ABC3.Found.Arakelov

open CategoryTheory

variable {A B C P : CommRingCat.{u}} {u : A ⟶ B} {v : A ⟶ C} {iB : B ⟶ P} {iC : C ⟶ P}

/-- ★★★★**押し出しで、片方の脚が全射なら反対の脚の核はイデアルの拡大である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `comap` のアフィン記述の核心である。 -/
theorem ker_of_isPushout (h : IsPushout u v iB iC) (hv : Function.Surjective v.hom) :
    RingHom.ker iB.hom = (RingHom.ker v.hom).map u.hom := by
  apply le_antisymm
  · set I : Ideal A := RingHom.ker v.hom with hIdef
    set J : Ideal B := I.map u.hom with hJdef
    let q : B ⟶ CommRingCat.of (B ⧸ J) := CommRingCat.ofHom (Ideal.Quotient.mk J)
    let e : C ≃+* (A ⧸ I) := (RingHom.quotientKerEquivOfSurjective hv).symm
    let lft : (A ⧸ I) →+* (B ⧸ J) :=
      Ideal.Quotient.lift I ((Ideal.Quotient.mk J).comp u.hom) (fun a ha =>
        Ideal.Quotient.eq_zero_iff_mem.2 (Ideal.mem_map_of_mem _ ha))
    let k : C ⟶ CommRingCat.of (B ⧸ J) := CommRingCat.ofHom (lft.comp e.toRingHom)
    have hcomm : u ≫ q = v ≫ k := by
      ext a
      show (Ideal.Quotient.mk J) (u.hom a) = lft (e (v.hom a))
      have he : e (v.hom a) = Ideal.Quotient.mk I a := by
        show (RingHom.quotientKerEquivOfSurjective hv).symm (v.hom a) = _
        rw [RingEquiv.symm_apply_eq]
        rfl
      rw [he]
      rfl
    have hd : iB ≫ h.desc q k hcomm = q := h.inl_desc _ _ _
    intro b hb
    have hq : q.hom b = 0 := by
      have := congrArg (fun (m : B ⟶ CommRingCat.of (B ⧸ J)) => m.hom b) hd
      rw [← this]
      show (h.desc q k hcomm).hom (iB.hom b) = 0
      rw [RingHom.mem_ker.1 hb, map_zero]
    exact Ideal.Quotient.eq_zero_iff_mem.1 hq
  · rw [Ideal.map_le_iff_le_comap]
    intro a ha
    have hw : iB.hom (u.hom a) = iC.hom (v.hom a) :=
      congrArg (fun (m : A ⟶ P) => m.hom a) h.w
    simp only [Ideal.mem_comap, RingHom.mem_ker] at ha ⊢
    rw [hw, ha, map_zero]


/-! ## ★出典の紐付け(`.src`) -/

def ker_of_isPushout.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——押し出しの核はイデアルの拡大)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
