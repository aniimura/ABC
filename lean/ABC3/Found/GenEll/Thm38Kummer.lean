/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Thm38Alpha
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★`q` が `K(ζ_l)` で `l` 乗にならないこと（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.20。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★これは何か

`§9-993` で `Theorem 3.8` の「`α` が像に入る」段の**行列表示**が取れた。
残るのは **Kummer の存在**——`l ∤ v_K(q)` なら
`σ(ζ) = ζ` かつ `σ(π) = ζ·π` なる `σ` があること——である。

Kummer 理論は「`Gal(K(ζ_l, q^{1/l})/K(ζ_l))` が `μ_l` と双対」なので、
**`σ` の存在は `q ∉ K(ζ_l)^{×l}` と同値**である。

★★★**本ファイルはその `q ∉ K(ζ_l)^{×l}` の側（付値の障害）を取る。**

## ★機構 —— 分岐指数が `l` と素であること

`K(ζ_l)/K` の分岐指数 `e` は `[K(ζ_l):K] ∣ l−1` を割るので **`l` と素**である
（`coprime_of_dvd_sub_one`）。したがって `q = y^l` なら

    `e · v_K(q) = l · v'(y)`

となり、`l` が素で `l ∤ e` だから **`l ∣ v_K(q)`**（`dvd_of_ramification_coprime`）。
★対偶が `not_lth_power_of_not_dvd_val` である。

## ★★残るのは Kummer 双対そのもの

☆`Gal(K(ζ_l, π)/K(ζ_l)) ↪ μ_l`（`σ ↦ σ(π)/π`）が**全射**であること
（`q ∉ K(ζ_l)^{×l}` のとき）。★これは体論の標準的な内容だが mathlib に
`IsCyclic`・`Polynomial.X_pow_sub_C` 周りの部品があるだけで、
Kummer 対応そのものは無い（2026-08-29 実測）。

★`.src` は条つき——指標には数えない。
-/

namespace ABC3.Found.GenEll

/-! ## ★★★★★分岐指数は `l` と素 -/

/-- ★★★★★**`e ∣ l−1` なら `e` は `l` と素**。

★`K(ζ_l)/K` の次数（したがって分岐指数）は `l−1` を割るので、`l` と素である。 -/
theorem coprime_of_dvd_sub_one {l e : ℕ} (hl : Nat.Prime l) (hepos : 0 < e)
    (he : e ∣ l - 1) : Nat.Coprime l e := by
  rw [Nat.Prime.coprime_iff_not_dvd hl]
  intro hdvd
  have hle : l ≤ e := Nat.le_of_dvd hepos hdvd
  have h2 : 1 < l := hl.one_lt
  have hsub : 0 < l - 1 := by omega
  have hle2 : e ≤ l - 1 := Nat.le_of_dvd hsub he
  omega

/-! ## ★★★★★★★付値の障害 -/

/-- ★★★★★★★**`e·v(q) = l·w` かつ `gcd(l,e) = 1` なら `l ∣ v(q)`**。

★`q = y^l` を `K(ζ_l)` の付値で読んだ形である（`e` は分岐指数）。 -/
theorem dvd_of_ramification_coprime {l e : ℕ} (hl : Nat.Prime l) (hcop : Nat.Coprime l e)
    {vq w : ℤ} (h : (e : ℤ) * vq = (l : ℤ) * w) : (l : ℤ) ∣ vq := by
  have hp : Prime (l : ℤ) := Nat.prime_iff_prime_int.mp hl
  have hdvd : (l : ℤ) ∣ (e : ℤ) * vq := ⟨w, h⟩
  rcases hp.dvd_mul.mp hdvd with hde | hdv
  · exfalso
    have hne : (l : ℕ) ∣ e := Int.natCast_dvd_natCast.mp hde
    have h1 : l ∣ Nat.gcd l e := Nat.dvd_gcd dvd_rfl hne
    rw [hcop] at h1
    exact hl.one_lt.ne' (Nat.dvd_one.mp h1)
  · exact hdv

/-- ★★★★★★★★★★★**`l ∤ v_K(q)` なら `q` は `K(ζ_l)` で `l` 乗にならない**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★付値で読めば「`e·v_K(q) = l·w` なる `w` は無い」である
（`e` は `K(ζ_l)/K` の分岐指数、`e ∣ l−1`）。
★★これが Kummer の存在（`σ(π) = ζ·π` なる `σ` がある）の**前提条件**である。 -/
theorem not_lth_power_of_not_dvd_val {l e : ℕ} (hl : Nat.Prime l) (hepos : 0 < e)
    (he : e ∣ l - 1) {vq : ℤ} (hnd : ¬ ((l : ℤ) ∣ vq)) (w : ℤ) :
    (e : ℤ) * vq ≠ (l : ℤ) * w := by
  intro h
  exact hnd (dvd_of_ramification_coprime hl (coprime_of_dvd_sub_one hl hepos he) h)

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def coprime_of_dvd_sub_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(K(ζ_l)/K の分岐指数は l と素)",
    sectionId := "genell-thm-3-8" }

def dvd_of_ramification_coprime.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(e·v(q) = l·w かつ gcd(l,e)=1 なら l ∣ v(q))",
    sectionId := "genell-thm-3-8" }

def not_lth_power_of_not_dvd_val.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l ∤ v_K(q) なら q は K(ζ_l) で l 乗にならない)",
    sectionId := "genell-thm-3-8" }

def not_lth_power_of_not_dvd_val.needs : List ABC3.Meta.ProofObligation :=
  [ .folklore
      ("Kummer 双対: `q ∉ K(ζ_l)^{×l}` なら " ++
       "Gal(K(ζ_l, q^{1/l})/K(ζ_l)) → μ_l (σ ↦ σ(π)/π) が**全射**であること。" ++
       "★体論の標準的な内容だが mathlib に Kummer 対応そのものは無い(2026-08-29 実測)。**残る**") 5,
    .folklore
      "[K(ζ_l):K] ∣ l−1(円分拡大の次数)——★これが分岐指数 e ∣ l−1 を与える" 3,
    .implicitStep
      ("★★★★★測定(2026-08-29): §9-993 で「α が像に入る」段の**行列表示**が取れ、" ++
       "本ファイルでその**前提条件(q ∉ K(ζ_l)^{×l})の付値の障害**が取れた。" ++
       "★機構は『K(ζ_l)/K の分岐指数 e は l−1 を割るので l と素』" ++
       "——q = y^l なら e·v_K(q) = l·v'(y) となり l ∣ v_K(q) になってしまう。" ++
       "★★残るのは **Kummer 双対そのもの**だけである") 6,
    .implicitStep
      ("★★Theorem 3.8 の残高(2026-08-29): " ++
       "(1) Kummer 双対(本ファイルで前提条件は取れた)、" ++
       "(2) l-巡回 ⟷ 安定直線の対応、(3) Lemma 3.7(Prop 3.4 待ち)、(4) torsionExt。" ++
       "★群論(Lemma 3.1, (iv)・§9-992・§9-993)と Galois 表現の構成(galRep)は済んでいる") 7 ]

end ABC3.Found.GenEll
