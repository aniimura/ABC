import ABC3.Found.GaloisRep.TowerSurj

/-!
# Galois (G5) 第 209 ブロック —— **★★★★★★★★逆極限の全射性と基底の生成性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★基底性の半分が閉じた

第 208 で塔の各段 `[l] : E[l^{m+1}] → E[l^m]` の全射が取れた。★本ブロックはそれを
**逆極限まで持ち上げる**——`tateProj n : T_l E → E[l^n]` が全射である。

### ★★★★★★両立列の構成

`R ∈ E[l^n]` から `l·f(m+1) = f(m)` を満たす列 `f` を作る。

| 段 | `m` | 定義 |
|---|---|---|
| 下 | `m ≤ n` | `f m := l^{n−m}·R`(`R` を割り戻す) |
| 上 | `m > n` | 各段の全射性で**再帰的に選ぶ**(`towerAux`) |

★上の段は `Exists.choose` の再帰であり、`towerAux hstep n R hR k` が
`E[l^{n+k}]` の元を返す。★★`towerSeq` はこの 2 つを `if m ≤ n` で貼り合わせたもので、
両立性 `l·f(m+1) = f(m)` は下・上・境界の 3 通りで確かめる。

### ★★★★★★★生成性

`e : T_l E ≃+ ℤ_l²` を通すと `tateVec n : ℤ_l² →+ E[l^n]` も全射になる。★任意の
`v : ℤ_l²` を `v = v₀·single₀ + v₁·single₁` と分解し、第 205 の
`φ(α·v) = (toZModPow n α).val · φ(v)` を当てると

    R = a·P + b·Q,   P := tateVec n single₀,  Q := tateVec n single₁

★★これが `hgen` の要求する**基底性のうち「生成する」半分**である。

## ★★残るのは `P` の位数が `l^n` ちょうどであることだけ

`hgen`(第 207)は「`e_{l^n}(P,Q)` が原始 `l^n` 乗根」——★非退化性(残件 (i)、
`normEDS` が楕円列であること、上流案件)と、`P, Q` が**基底**であること。
★★本ブロックで生成性が済んだので、残るのは自由性(位数がちょうど `l^n`)である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `towerAux` ほか | 塔を上に伸ばす再帰 |
| `towerSeq_compat` | ★★★★★★両立列であること |
| `tateProj_surjective` | ★★★★★★★**`T_l E → E[l^n]` は全射** |
| `addHom_repr` | ★★★★★★任意のベクトルを基底で読む |
| `exists_smul_repr` | ★★★★★★★**`P, Q` は `E[l^n]` を生成する** |
-/

namespace ABC3.Found.GaloisRep

open ABC3.Interface.GaloisRep WeierstrassCurve WeierstrassCurve.Affine

universe u

variable {A : Type u} [AddCommGroup A] {l : ℕ}

/-- 塔を上に伸ばす再帰——`towerAux hstep n R hR k ∈ A[l^{n+k}]` である。 -/
noncomputable def towerAux
    (hstep : ∀ (m : ℕ) (T : A), (l ^ m) • T = 0 → ∃ S : A, (l ^ (m + 1)) • S = 0 ∧ l • S = T)
    (n : ℕ) (R : A) (hR : (l ^ n) • R = 0) :
    (k : ℕ) → {P : A // (l ^ (n + k)) • P = 0}
  | 0 => ⟨R, by rw [Nat.add_zero]; exact hR⟩
  | (k + 1) =>
    ⟨(hstep (n + k) (towerAux hstep n R hR k).1 (towerAux hstep n R hR k).2).choose, by
      have h := (hstep (n + k) (towerAux hstep n R hR k).1
        (towerAux hstep n R hR k).2).choose_spec.1
      rw [show n + (k + 1) = (n + k) + 1 by omega]
      exact h⟩

theorem towerAux_zero
    (hstep : ∀ (m : ℕ) (T : A), (l ^ m) • T = 0 → ∃ S : A, (l ^ (m + 1)) • S = 0 ∧ l • S = T)
    (n : ℕ) (R : A) (hR : (l ^ n) • R = 0) :
    (towerAux hstep n R hR 0).1 = R := rfl

theorem towerAux_succ
    (hstep : ∀ (m : ℕ) (T : A), (l ^ m) • T = 0 → ∃ S : A, (l ^ (m + 1)) • S = 0 ∧ l • S = T)
    (n : ℕ) (R : A) (hR : (l ^ n) • R = 0) (k : ℕ) :
    l • (towerAux hstep n R hR (k + 1)).1 = (towerAux hstep n R hR k).1 :=
  (hstep (n + k) (towerAux hstep n R hR k).1 (towerAux hstep n R hR k).2).choose_spec.2

/-- 両立列——`m ≤ n` では `R` を割り戻し、`m > n` では塔を伸ばす。 -/
noncomputable def towerSeq
    (hstep : ∀ (m : ℕ) (T : A), (l ^ m) • T = 0 → ∃ S : A, (l ^ (m + 1)) • S = 0 ∧ l • S = T)
    (n : ℕ) (R : A) (hR : (l ^ n) • R = 0) (m : ℕ) : A :=
  if m ≤ n then (l ^ (n - m)) • R else (towerAux hstep n R hR (m - n)).1

theorem towerSeq_torsion
    (hstep : ∀ (m : ℕ) (T : A), (l ^ m) • T = 0 → ∃ S : A, (l ^ (m + 1)) • S = 0 ∧ l • S = T)
    (n : ℕ) (R : A) (hR : (l ^ n) • R = 0) (m : ℕ) :
    (l ^ m) • towerSeq hstep n R hR m = 0 := by
  rw [towerSeq]
  split_ifs with hm
  · rw [← mul_smul, ← pow_add, show m + (n - m) = n by omega, hR]
  · have h := (towerAux hstep n R hR (m - n)).2
    rw [show l ^ m = l ^ (n + (m - n)) by rw [show n + (m - n) = m by omega]]
    exact h

theorem towerSeq_at
    (hstep : ∀ (m : ℕ) (T : A), (l ^ m) • T = 0 → ∃ S : A, (l ^ (m + 1)) • S = 0 ∧ l • S = T)
    (n : ℕ) (R : A) (hR : (l ^ n) • R = 0) :
    towerSeq hstep n R hR n = R := by
  rw [towerSeq, if_pos le_rfl, Nat.sub_self, pow_zero, one_smul]

/-- ★★★★★★**両立列であること**——`l·f(m+1) = f(m)`。 -/
theorem towerSeq_compat
    (hstep : ∀ (m : ℕ) (T : A), (l ^ m) • T = 0 → ∃ S : A, (l ^ (m + 1)) • S = 0 ∧ l • S = T)
    (n : ℕ) (R : A) (hR : (l ^ n) • R = 0) (m : ℕ) :
    l • towerSeq hstep n R hR (m + 1) = towerSeq hstep n R hR m := by
  rw [towerSeq, towerSeq]
  by_cases hm : m + 1 ≤ n
  · rw [if_pos hm, if_pos (by omega : m ≤ n), ← mul_smul, ← pow_succ',
      show n - (m + 1) + 1 = n - m by omega]
  · rw [if_neg hm]
    by_cases hm2 : m ≤ n
    · rw [if_pos hm2, show n - m = 0 by omega, pow_zero, one_smul,
        show m + 1 - n = 0 + 1 by omega]
      have h := towerAux_succ hstep n R hR 0
      rw [h, towerAux_zero]
    · rw [if_neg hm2, show m + 1 - n = (m - n) + 1 by omega]
      exact towerAux_succ hstep n R hR (m - n)

/-- ★★★★★★★**`tateProj n : T_l E → E[l^n]` は全射である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★塔の各段の全射性(第 208)を逆極限まで持ち上げた形である。 -/
theorem tateProj_surjective {F : Type} [Field F] [DecidableEq F]
    (W : WeierstrassCurve F) (l n : ℕ)
    (hstep : ∀ (m : ℕ) (T : W.toAffine.Point), (l ^ m) • T = 0 →
      ∃ S : W.toAffine.Point, (l ^ (m + 1)) • S = 0 ∧ l • S = T)
    (R : torsionPoints W (l ^ n)) :
    ∃ f : tateModule W l, tateProj W l n f = R := by
  refine ⟨⟨fun m => ⟨towerSeq hstep n (R : W.toAffine.Point) R.2 m,
      towerSeq_torsion hstep n (R : W.toAffine.Point) R.2 m⟩,
    fun m => towerSeq_compat hstep n (R : W.toAffine.Point) R.2 m⟩, ?_⟩
  exact Subtype.ext (towerSeq_at hstep n (R : W.toAffine.Point) R.2)

/-- ★★★★★★**任意のベクトルを基底で読む**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 206 の `addHom_mulVec_single` の、行列を経由しない版である。 -/
theorem addHom_repr {A : Type u} [AddCommGroup A] {l : ℕ} [Fact l.Prime]
    (n : ℕ) (φ : (Fin 2 → ℤ_[l]) →+ A) (htor : ∀ v, (l ^ n) • φ v = 0)
    (v : Fin 2 → ℤ_[l]) :
    φ v = ((PadicInt.toZModPow n (v 0)).val) • φ (Pi.single 0 1)
      + ((PadicInt.toZModPow n (v 1)).val) • φ (Pi.single 1 1) := by
  have hdec : v = v 0 • (Pi.single 0 1 : Fin 2 → ℤ_[l])
      + v 1 • (Pi.single 1 1 : Fin 2 → ℤ_[l]) := by
    funext i
    fin_cases i <;> simp
  conv_lhs => rw [hdec]
  rw [map_add, addHom_padic_smul_of_torsion n φ htor,
    addHom_padic_smul_of_torsion n φ htor]

/-- ★★★★★★★**`P, Q` は `E[l^n]` を生成する**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 207 の `hgen` が要求する**基底性のうち「生成する」半分**である。 -/
theorem exists_smul_repr {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L]
    [Algebra K L] (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (n : ℕ)
    (hstep : ∀ (m : ℕ) (T : (W.baseChange L).toAffine.Point), (l ^ m) • T = 0 →
      ∃ S : (W.baseChange L).toAffine.Point, (l ^ (m + 1)) • S = 0 ∧ l • S = T)
    {R : (W.baseChange L).toAffine.Point} (hR : (l ^ n) • R = 0) :
    ∃ a b : ℕ, R = a • tateVec W l e n (Pi.single 0 1)
      + b • tateVec W l e n (Pi.single 1 1) := by
  obtain ⟨f, hf⟩ := tateProj_surjective (W.baseChange L) l n hstep ⟨R, hR⟩
  have hv : tateVec W l e n (e f) = R := by
    rw [tateVec_apply, AddEquiv.symm_apply_apply, hf]
  refine ⟨((PadicInt.toZModPow n ((e f) 0)).val), ((PadicInt.toZModPow n ((e f) 1)).val), ?_⟩
  rw [← hv]
  exact addHom_repr n (tateVec W l e n) (tateVec_torsion W l e n) (e f)

/-! ## ★出典の紐付け(`.src`) -/

def tateProj_surjective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 進表現の行列式——Tate 加群から E[l^n] への射影の全射性)",
    sectionId := "genell-thm-3-8" }

def addHom_repr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 進表現の行列式——ベクトルを基底で読むこと)",
    sectionId := "genell-thm-3-8" }

def exists_smul_repr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l 進表現の行列式——P, Q が E[l^n] を生成すること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
