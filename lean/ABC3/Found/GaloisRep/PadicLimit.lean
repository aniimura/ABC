import ABC3.Found.GaloisRep.IndepTower
import Mathlib.NumberTheory.Padics.RingHoms

/-!
# Galois (G2) 第 75 ブロック —— **★★★`ℤ_p` は `ZMod (p^n)` の逆極限**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★mathlib にある半分と、無い半分

mathlib は `PadicInt.toZModPow : ℤ_p →+* ZMod (p^n)` と、その**単射性**
(`ext_of_toZModPow`)を持つ。★**全射性**——両立する剰余の列が `ℤ_p` から来ること——は
直接の形では無い。

★★しかし `PadicInt.lift` がある:任意の環 `R` から `ZMod (p^n)` への
両立する射の族は `R →+* ℤ_p` を誘導する。

★★★そこで **`R` として「両立する列のなす部分環」そのもの**を取れば、
`lift` が求める全射性を与える——これが本ブロックの仕掛けである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `compatSubring` | ★両立する剰余の列のなす部分環 |
| `exists_padicInt_of_compat` | ★★★**逆極限表示**(全射性) |
-/

namespace ABC3.Found.GaloisRep

/-- ★★**両立する剰余の列**のなす部分環。 -/
def compatSubring (p : ℕ) : Subring (∀ n : ℕ, ZMod (p ^ n)) where
  carrier := {x | ∀ (m n : ℕ) (h : m ≤ n),
    ZMod.castHom (pow_dvd_pow p h) (ZMod (p ^ m)) (x n) = x m}
  mul_mem' := by
    intro a b ha hb m n h
    simp only [Pi.mul_apply, map_mul, ha m n h, hb m n h]
  one_mem' := by intro m n h; exact map_one _
  add_mem' := by
    intro a b ha hb m n h
    simp only [Pi.add_apply, map_add, ha m n h, hb m n h]
  zero_mem' := by intro m n h; exact map_zero _
  neg_mem' := by
    intro a ha m n h
    simp only [Pi.neg_apply, map_neg, ha m n h]

/-- ★★★**`ℤ_p` は `ZMod (p^n)` の逆極限である**——両立する列は `p` 進整数から来る。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem exists_padicInt_of_compat (p : ℕ) [Fact p.Prime] (x : ∀ n : ℕ, ZMod (p ^ n))
    (hx : ∀ n : ℕ, ZMod.castHom (pow_dvd_pow p (Nat.le_succ n)) (ZMod (p ^ n)) (x (n + 1)) = x n) :
    ∃ α : ℤ_[p], ∀ n, PadicInt.toZModPow n α = x n := by
  have hfull : ∀ (m n : ℕ) (h : m ≤ n),
      ZMod.castHom (pow_dvd_pow p h) (ZMod (p ^ m)) (x n) = x m := by
    intro m n h
    induction n, h using Nat.le_induction with
    | base =>
      have hid : ZMod.castHom (pow_dvd_pow p (le_refl m)) (ZMod (p ^ m)) = RingHom.id _ :=
        RingHom.ext_zmod _ _
      rw [hid]; rfl
    | succ n hmn ih =>
      have hcomp : ZMod.castHom (pow_dvd_pow p (le_trans hmn (Nat.le_succ n))) (ZMod (p ^ m))
          = (ZMod.castHom (pow_dvd_pow p hmn) (ZMod (p ^ m))).comp
            (ZMod.castHom (pow_dvd_pow p (Nat.le_succ n)) (ZMod (p ^ n))) :=
        RingHom.ext_zmod _ _
      rw [hcomp]
      simp only [RingHom.coe_comp, Function.comp_apply, hx n, ih]
  let f : (k : ℕ) → (compatSubring p) →+* ZMod (p ^ k) := fun k =>
    (Pi.evalRingHom _ k).comp ((compatSubring p).subtype)
  have f_compat : ∀ (k1 k2 : ℕ) (hk : k1 ≤ k2),
      (ZMod.castHom (pow_dvd_pow p hk) (ZMod (p ^ k1))).comp (f k2) = f k1 := by
    intro k1 k2 hk
    ext y
    exact y.2 k1 k2 hk
  refine ⟨PadicInt.lift f_compat ⟨x, hfull⟩, fun n => ?_⟩
  exact RingHom.congr_fun (PadicInt.lift_spec f_compat n) ⟨x, hfull⟩

/-! ## ★出典の紐付け(`.src`) -/

def exists_padicInt_of_compat.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Tate 加群——Z_p の逆極限表示)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
