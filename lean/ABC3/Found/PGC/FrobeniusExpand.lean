import Mathlib.RingTheory.MvPowerSeries.Expand
import Mathlib.FieldTheory.Finite.Basic

/-!
# 有限体上の冪級数の Frobenius 恒等式(`sorry` 無し)

Lubin-Tate 理論(局所類体論の構成、`ResearchPaper/blocked-leaves.json` の
`[pGC] Proposition 1.2 / ...`)の**存在補題**(`f, g ∈ 𝔉_π` に対し `F(X,Y)` を
次数ごとに再帰的に構成する)の核心にある可除性(「次数 n+1 の剰余項 `R_n` が
`π` で割り切れる」)は、剰余体 `𝓀 = O/π`(位数 `q = p^f` の有限体)上で

  `φ(X)^q ≡ φ(X^q) (mod π)`(`φ` は `𝓀` 係数の冪級数)

という Frobenius 型の恒等式に帰着する(§9-1467 台で確認・記録した議論、
`memory/padic-log-additivity-blocked.md` 参照)。**本ファイルはこの恒等式を
一般の有限体について sorry 無しで証明した。**

## 鍵になった mathlib の道具(実測、2026-09-04)

- `MvPowerSeries.expand (k) (hk : k≠0) : MvPowerSeries σ R →ₐ[R] MvPowerSeries σ R`
  ——`Xᵢ ↦ Xᵢᵏ` の代入(`expand_X` で確認)。
- `MvPowerSeries.map_iterateFrobenius_expand (p) (hp) (f) (n) :
     map (iterateFrobenius R p n) (expand (p^n) _ f) = f ^ p^n`
  ——`RingTheory/MvPowerSeries/Expand.lean` に**既にある**。探すまで存在を知らなかった
  (`grep` では見つからず、`frobenius.*subst`/`subst.*frobenius` 系のキーワードでも
  当たらなかった——`expand`という語を知らないと到達できない命名だった)。
- 有限体 `F`(位数 `q = p^f`)では `iterateFrobenius F p f = id`(`FiniteField.pow_card`
  ——`∀ x:F, x^(Fintype.card F) = x` から)。この2つを合わせると
  `φ^q = map(id)(expand q φ) = expand q φ` が**そのまま出る**。

## まだ無いもの

本補題は Lubin-Tate 構成の可除性議論の**核**だが、そこから実際に「次数ごとに
係数を再帰的に決定する」構成全体(整合する係数列の定義、無限次までの検証)は
別途大きな作業として残る——`ResearchPaper/blocked-leaves.json` の
`[pGC] Proposition 1.2 / ...` を参照。
-/

namespace ABC3.Found.PGC

/-- **有限体上の冪級数の Frobenius 恒等式**: 位数 `q = p^f` の有限体 `F` 上の
任意の多変数冪級数 `φ` について、`φ` を `q` 乗したものは、各変数を `q` 乗で
置き換えたもの(`MvPowerSeries.expand`)に一致する。

Lubin-Tate 構成で `f, g ∈ 𝔉_π` の可除性議論の核心(`F(X)^q ≡ F(X^q) (mod π)`)
はこの恒等式を `F := O/π`(有限体)に適用したもの。 -/
theorem mvPowerSeries_pow_card_eq_expand {σ F : Type*} [Field F] [Fintype F] {p : ℕ}
    [ExpChar F p] {f : ℕ} (hq : Fintype.card F = p ^ f) (φ : MvPowerSeries σ F) :
    φ ^ (p ^ f) = MvPowerSeries.expand (p ^ f) (pow_ne_zero f (expChar_ne_zero F p)) φ := by
  have hid : iterateFrobenius F p f = RingHom.id F := by
    ext x
    show x ^ p ^ f = x
    rw [← hq]
    exact FiniteField.pow_card x
  have h := MvPowerSeries.map_iterateFrobenius_expand p (expChar_ne_zero F p) φ f
  rw [hid] at h
  simpa using h.symm

end ABC3.Found.PGC
