import ABC3.Found.Arakelov.Def11Complete

/-!
# 新しい `deg_F` の**有限素点側は退化していない** —— `Definition 1.2` を上げるための検査

`Check/GenEll/Def12Degenerate.lean` は、**古い**設計（`Found/Arakelov/ADegEmb.lean`）の

    degFOf F x = -(∑_{σ : F →+* ℂ} x.2.fn (embPoint F σ)) / [F:ℚ]

が `x.1`（`Pic` 類）を式に含まないこと（`ht_indep_pic`、`rfl` で落ちる）を機械にかけ、
それを理由に `Definition 1.2` の項目全体の `.src` を下げた（2026-08-27）。

★**本ファイルはその逆向きの検査である**——`§9-776`〜`§9-782` で作った**新しい**

    degArithPre L s = degFinPre L s + archDeg L s
    degFinPre L s   = log #(Γ_pre(L) / Γ(X,⊤)·s)

の**有限素点側が恒等的に 0 ではない**ことを示す。

**これは原典の主張ではない**（我々のモデルについての事実）ので `.src` を持たない。

## ★★★何を示すか

    degFinPre L (2·s) ≠ degFinPre L s        （`degFinPre_two_ne`）
    degFinPre L (2·s) ≠ 0 ∨ degFinPre L s ≠ 0（`exists_degFinPre_ne_zero`）

★古い設計の退化（有限素点側の寄与が**無い**）であれば、どちらも成り立たない。
★★したがってこれは古い失敗形の**機械による反証**である。

## ★★★★機構

`degFinPre_smul`（`§9-780`）:

    degFinPre L (c·s) = degFinPre L s + log #(Γ(X,⊤)/(c))

に `c = 2` を入れ、`#(𝓞_F/(2)) = |N(2)| = 2^{[F:ℚ]} > 1` を使うだけである。

★`Ideal.absNorm I = Nat.card (R ⧸ I)` は**定義的に等しい**ので、
`Ideal.absNorm_span_singleton` と `Algebra.norm_algebraMap` がそのまま効く。

## ★★★★★これで何が言えて、何が言えないか

★**言えること**: 新しい `deg_F` は有限素点側を実際に測っている
——古い `ht_indep_pic` 型の退化は**起きていない**。

★★**言えないこと**: 「`deg_F` が `APic` 上で単射である」等の強い非退化性。
それは類数や単数の情報を要し、本項目（`Definition 1.2` は `ht` の**命名**である）には要らない。
-/

namespace ABC3.Check.GenEll

open AlgebraicGeometry CategoryTheory Opposite NumberField
open ABC3.Found.Arakelov

/-! ## ★★★`𝓞_F/(2)` の位数 -/

/-- ★★`#(𝓞_F/(2)) = 2^{[F:ℚ]}`。

★`Ideal.absNorm I = Nat.card (R ⧸ I)` は定義的に等しい。 -/
theorem card_quotient_two_ringOfIntegers (F : Type) [Field F] [NumberField F] :
    Nat.card ((𝓞 F) ⧸ (Ideal.span {(2 : 𝓞 F)})) = 2 ^ (Module.finrank ℚ F) := by
  show Ideal.absNorm (Ideal.span {(2 : 𝓞 F)}) = _
  rw [Ideal.absNorm_span_singleton]
  have h2 : (2 : 𝓞 F) = algebraMap ℤ (𝓞 F) 2 := by norm_num
  rw [h2, Algebra.norm_algebraMap, RingOfIntegers.rank]
  simp

/-- ★★`Γ(Spec 𝓞_F, ⊤)` の側でも同じ。 -/
theorem card_quotient_two_gammaSpec (F : Type) [Field F] [NumberField F] :
    Nat.card ((Γ(Spec (CommRingCat.of (𝓞 F)),
        (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)
      ⧸ (Ideal.span {(2 : (Γ(Spec (CommRingCat.of (𝓞 F)),
          (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type))}))
      = 2 ^ (Module.finrank ℚ F) := by
  rw [card_quotient_gammaSpec (CommRingCat.of (𝓞 F)) 2]
  have h : (Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).hom.hom 2 = (2 : 𝓞 F) :=
    map_ofNat (Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).hom.hom 2
  rw [h]
  exact card_quotient_two_ringOfIntegers F

theorem two_ne_zero_gammaSpec (F : Type) [Field F] [NumberField F] :
    (2 : (Γ(Spec (CommRingCat.of (𝓞 F)),
      (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)) ≠ 0 := by
  intro h
  have h2 : (gammaSpecRingEquiv (CommRingCat.of (𝓞 F))) 2
      = (gammaSpecRingEquiv (CommRingCat.of (𝓞 F))) 0 := by rw [h]
  rw [map_ofNat, map_zero] at h2
  exact (two_ne_zero : (2 : 𝓞 F) ≠ 0) h2

/-! ## ★★★★★★★★有限素点側は退化していない -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★**新しい `deg_F` の有限素点側は切断を実際に見ている**。

    `degFinPre L (2·s) ≠ degFinPre L s`

★古い設計の退化（有限素点側の寄与が無い）ならこれは成り立たない。 -/
theorem degFinPre_two_ne (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F))))
    (s : (L.carrier.sheaf.obj (op ⊤) : Type)) (hs : s ≠ 0) :
    degFinPre L ((2 : (Γ(Spec (CommRingCat.of (𝓞 F)),
        (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)) • s) ≠ degFinPre L s := by
  rw [degFinPre_smul L s hs 2 (two_ne_zero_gammaSpec F)
      (fun r hr => finite_quotient_gammaSpec F r hr),
    card_quotient_two_gammaSpec F]
  have hn : 0 < Module.finrank ℚ F := Module.finrank_pos
  have h1 : (1 : ℝ) < ((2 ^ (Module.finrank ℚ F) : ℕ) : ℝ) := by
    have : (1 : ℕ) < 2 ^ (Module.finrank ℚ F) := Nat.one_lt_two_pow_iff.2 hn.ne'
    exact_mod_cast this
  have hlog : 0 < Real.log ((2 ^ (Module.finrank ℚ F) : ℕ) : ℝ) := Real.log_pos h1
  intro hcon
  exact absurd (by linarith [hcon] : Real.log ((2 ^ (Module.finrank ℚ F) : ℕ) : ℝ) = 0)
    hlog.ne'

/-- ★★★★★★★★★**有限素点側は恒等的に `0` ではない**
—— 古い `ht_indep_pic` 型の退化の**反証**。 -/
theorem exists_degFinPre_ne_zero (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F))))
    (s : (L.carrier.sheaf.obj (op ⊤) : Type)) (hs : s ≠ 0) :
    degFinPre L ((2 : (Γ(Spec (CommRingCat.of (𝓞 F)),
        (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)) • s) ≠ 0
      ∨ degFinPre L s ≠ 0 := by
  by_contra hcon
  push_neg at hcon
  exact degFinPre_two_ne F L s hs (hcon.1.trans hcon.2.symm)

end ABC3.Check.GenEll
