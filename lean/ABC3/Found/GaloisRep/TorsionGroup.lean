import ABC3.Found.GaloisRep.TorsionCard

/-!
# Galois (G1) 第 67 ブロック —— **★★★★★独立な 2 元から `(ℤ/n)²`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★`structure_eq` の設計を確定する

`Interface/GaloisRep/Torsion.lean` の `structure_eq` は現状

    Nonempty (torsionPoints W n ≃ (ZMod n × ZMod n))

と **`Equiv`(型の全単射)** で書いてある。★これは第 66 ブロックの `Nat.card = n²`
だけで埋まってしまう——★★**原典 Theorem 3.8 の主張(群同型)より弱い**。

★★★したがって `Interface` を **`≃+`(`AddEquiv`)** に強めることにする(§9-400 に記録)。
本ブロックはその「易しい半分」である。

## ★★易しい半分

`n•a = n•b = 0` で **独立**(`i•a + j•b = 0 ⟹ n ∣ i ∧ n ∣ j`)な 2 元があれば、

    (ℤ/n)² →+ A,  (i, j) ↦ i•a + j•b

は単射で、`#A = n²` から**全単射**になる。

## ★★残る半分(§9-400)

独立な 2 元の存在。★素数冪 `n = p^k` では数え上げで出る:

| 集合 | 個数 |
|---|---|
| 位数 `p^k` の元 | `p^{2k} − p^{2k−2}` |
| `⟨a⟩ ∩ ⟨b⟩ ≠ 0` となる `b` | `≤ p^{2k−1}` |

★★`p² − 1 > p` なので前者が大きい ✅
★★★一般の `n` は中国剰余定理で素数冪へ落とす。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `torsionSub` | ★`n` 捩れ部分群 |
| `addEquiv_of_indep` | ★★★★★**独立な 2 元 ⟹ `(ℤ/n)² ≃+ A`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- ★**`n` 捩れ部分群**。 -/
def torsionSub (W : WeierstrassCurve F) (n : ℕ) : AddSubgroup W.toAffine.Point where
  carrier := {P | n • P = 0}
  add_mem' := by intro a b ha hb; simp only [Set.mem_setOf_eq] at *; rw [smul_add, ha, hb, add_zero]
  zero_mem' := by simp
  neg_mem' := by intro a ha; simp only [Set.mem_setOf_eq] at *; rw [smul_neg, ha, neg_zero]

theorem mem_torsionSub {n : ℕ} {P : W.toAffine.Point} :
    P ∈ torsionSub W n ↔ n • P = 0 := Iff.rfl

theorem torsionSub_smul (n : ℕ) (a : torsionSub W n) : n • a = 0 := by
  ext
  exact a.2

/-- ★★★★★**独立な位数 `n` の 2 元があり `#A = n²` なら `(ℤ/n)² ≃+ A`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem addEquiv_of_indep {A : Type u} [AddCommGroup A] [Finite A] (n : ℕ) (hn : 1 ≤ n)
    (a b : A) (ha : n • a = 0) (hb : n • b = 0)
    (hind : ∀ i j : ℤ, i • a + j • b = 0 → ((n : ℤ) ∣ i ∧ (n : ℤ) ∣ j))
    (hcard : Nat.card A = n ^ 2) :
    Nonempty ((ZMod n × ZMod n) ≃+ A) := by
  haveI : NeZero n := ⟨by omega⟩
  have hza : (zmultiplesHom A a) (n : ℤ) = 0 := by simpa using ha
  have hzb : (zmultiplesHom A b) (n : ℤ) = 0 := by simpa using hb
  set fa : ZMod n →+ A := ZMod.lift n ⟨zmultiplesHom A a, hza⟩ with hfa
  set fb : ZMod n →+ A := ZMod.lift n ⟨zmultiplesHom A b, hzb⟩ with hfb
  have hfaval : ∀ i : ZMod n, fa i = (i.val : ℤ) • a := by
    intro i
    have hc : ((i.val : ℤ) : ZMod n) = i := by simp
    rw [hfa, ← hc, ZMod.lift_coe]
    simp
  have hfbval : ∀ i : ZMod n, fb i = (i.val : ℤ) • b := by
    intro i
    have hc : ((i.val : ℤ) : ZMod n) = i := by simp
    rw [hfb, ← hc, ZMod.lift_coe]
    simp
  set f : (ZMod n × ZMod n) →+ A :=
    (fa.comp (AddMonoidHom.fst (ZMod n) (ZMod n))) + (fb.comp (AddMonoidHom.snd (ZMod n) (ZMod n)))
      with hf
  have hfval : ∀ p : ZMod n × ZMod n, f p = (p.1.val : ℤ) • a + (p.2.val : ℤ) • b := by
    intro p; rw [hf]; simp [hfaval, hfbval]
  have hinj : Function.Injective f := by
    rw [injective_iff_map_eq_zero]
    rintro ⟨i, j⟩ hij
    rw [hfval] at hij
    obtain ⟨hd1, hd2⟩ := hind _ _ hij
    refine Prod.ext ?_ ?_ <;> simp only
    · have hz : ((i.val : ℤ) : ZMod n) = 0 := by
        rw [ZMod.intCast_zmod_eq_zero_iff_dvd]; exact_mod_cast hd1
      simpa using hz
    · have hz : ((j.val : ℤ) : ZMod n) = 0 := by
        rw [ZMod.intCast_zmod_eq_zero_iff_dvd]; exact_mod_cast hd2
      simpa using hz
  have hbij : Function.Bijective f := by
    rw [Nat.bijective_iff_injective_and_card]
    refine ⟨hinj, ?_⟩
    rw [hcard, Nat.card_prod, Nat.card_zmod]
    ring
  exact ⟨AddEquiv.ofBijective f hbij⟩

/-! ## ★出典の紐付け(`.src`) -/

def addEquiv_of_indep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(独立な 2 元から (Z/n)² への同型)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
