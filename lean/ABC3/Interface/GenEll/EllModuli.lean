import ABC3.Interface.GenEll.GaloisRep
import Mathlib.Data.Real.Basic
import Mathlib.Data.Set.Finite.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic

/-!
# [GenEll] §3–§4 大域理論 —— `M_ell(ℚ̄)` 上の高さの `Interface`

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.17–p.23。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★`TorsionGaloisRepData` を**拡張する**(並立させない)

★`Interface/GenEll/GaloisRep.lean` の `TorsionGaloisRepData` は
`Theorem 3.8` の statement のためだけに作った**狭い**構造体である。
§3 の残りと §4 は同じ語(`M_ell(ℚ̄)`・`ht^Falt`・compactly bounded・Galois-finite)を
使うので、**別に posit すると 2 つの語彙が並立して腐る**。

★ゆえに `extends` で**機械的に接続する**。
`Theorem 3.8` の statement と `Corollary 4.3/4.4` の statement は、
**同じ `EllClass`・同じ `faltingsHeight` の上に立つことが型で保証される**。

## ★mathlib 実測(2026-08-16)

| 語 | mathlib |
|---|---|
| Faltings 高さ `ht^Falt` | ★**0 件** |
| `deg_∞` / `ht_∞`(無限遠因子) | ★**0 件** |
| `log-diff_{M̄_ell}` | ★**0 件**(`Definition 1.5, (iii)` を `M̄_ell` に適用したもの) |
| 楕円曲線のモジュライスタック `M̄_ell` | ★**0 件** |

## ★★★★★★★★ 2026-08-26 の訂正(欠陥 #8 を塞ぐ)

以前は本構造体は**公理を 1 つも持たなかった**。
★そのため `Skeleton/GenEll/Section3.lean` の `prop_3_4` は
`deg_∞ ≡ 0`・`ht_∞ = n`・`ht^Falt ≡ 0` などの退化した `D` で**破れていた**
(`Check/GenEll/Section3NotProvable.lean` で実測)。

★★そこで原文が `Proposition 3.4` の証明で**実際に引いている 4 つ**を欄に出した:

| 欄 | 原文の典拠 |
|---|---|
| `degInf_le_htInf` | [GenEll] `Proposition 1.6` の証明 |
| `htInf_bdeq_faltings` | [Silv2] `Proposition 2.1` + [FC] Ch. V, `Proposition 4.5` |
| `faltingsHeight_bddBelow` | 古典的(Faltings 高さは下に有界) |
| `northcott` | [GenEll] `Proposition 1.4, (iv)` |

★★★**`ε` の入った 2 つの不等式はこれらから導かれる**——
posit しているのは原文が引く外部文献と内部参照だけである。
-/

namespace ABC3.Interface.GenEll

open ABC3.Meta

/-- **[GenEll] §3–§4 の大域理論**を受ける `Interface`。

★`TorsionGaloisRepData`(`Theorem 3.8` 用)を**拡張**している。
追加したフィールドはすべて原文 p.17–p.23 に**実際に現れる語**である。 -/
structure EllModuliData extends TorsionGaloisRepData where
  /-- 原文 `deg_∞`(無限遠因子の次数、`Lemma 3.2, (ii)` の大域版)。 -/
  degInf : EllClass → ℝ
  /-- 原文 `ht_∞`(無限遠因子に付随する高さ)。 -/
  htInf : EllClass → ℝ
  /-- 原文 `log-diff_{M̄_ell}`(`Definition 1.5, (iii)` を `M̄_ell` に適用したもの)。 -/
  logDiffMell : EllClass → ℝ
  /-- 原文 `M_ell(ℚ̄)^{≤d}`(`Example 1.3, (i)`)。 -/
  degLe : ℕ → Set EllClass
  /-- 原文「`E_L` … with **semi-stable reduction** at all the finite primes of `L`」。 -/
  SemiStable : Curve → Prop
  /-- 原文「`E_L` admits an **l-cyclic** subgroup scheme `H_L ⊆ E_L`」。 -/
  HasLCyclic : Curve → ℕ → Prop
  /-- 原文「`L` is a **minimal field of definition** of the point `[E_L]`」。 -/
  MinimalField : Curve → Prop
  /-- 原文「The Galois representation `Gal(ℚ̄/L) → GL₂(ℤ_l)` associated to `E_L`
  is **surjective**」(`Corollary 4.3/4.4`, (b))。 -/
  ImageSurjective : Curve → ℕ → Prop
  /-- 原文「`l_•` is **prime** to the primes of `ℚ` that **ramify** in `L`, as well as to
  the **ramification indices** of primes of `ℚ` in `L`」(`Corollary 4.3/4.4`, (a))。 -/
  PrimeToRamification : Curve → ℕ → Prop
  /-- 原文「`E_L` has at least **one** prime of [bad] **multiplicative** reduction」。
  ★`TorsionGaloisRepData.HasPotMultRed`(**潜在的**乗法還元)とは**別の述語**である
  ——`Lemma 3.7` は両者を使い分けている(p.18 目視確認)。 -/
  HasMultRed : Curve → Prop
  /-- 原文「`l_∘`, `l_•` are **prime** to the **primes of potentially multiplicative
  reduction**」(`Corollary 4.3/4.4`, (a) の前半)。
  ★`PrimeToLocalHeights`(局所高さと素)とは**別の条件**である——
  原文 (a) は "as well as to" で 2 つを並べている(p.22 目視確認)。 -/
  PrimeToMultPrimes : Curve → ℕ → Prop
  /-- ★★★★★★**`deg_∞ ≲ ht_∞`**——`Proposition 3.4` の最初の `≲`。

  原文 (GenEll p.17):
  > Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

  ★原文はこれを `Proposition 1.6` の証明から取る。
  ★★**向きは `deg_∞(x) − ht_∞(x) ≤ C`**(abc の向き、`BDge`)である——
  逆向きにすると `ht_∞` が上に有界になる(`Check/GenEll/Section3NotProvable.lean`)。 -/
  degInf_le_htInf : ∃ C : ℝ, ∀ x, degInf x - htInf x ≤ C
  /-- ★★★★★★★**`ht_∞ ≈ 12·ht^Falt`**——無限遠での計量の対数的特異性。

  原文 (GenEll p.17):
  > Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

  ★原文の典拠は **[Silv2] Proposition 2.1** と **[FC] Ch. V, Proposition 4.5**。
  ★★★これが `Proposition 3.4` の 2 番目・3 番目の `≲` を**導く**。 -/
  htInf_bdeq_faltings : ∃ C : ℝ, ∀ x, |htInf x - 12 * faltingsHeight x| ≤ C
  /-- ★★★★**Faltings 高さは下に有界**。

  ★`ε` の入った不等式を出すのに要る——`12ε·ht^Falt` を下から抑えねばならない。 -/
  faltingsHeight_bddBelow : ∃ B : ℝ, ∀ x, B ≤ faltingsHeight x
  /-- ★★★★★★**Northcott の有限性**——原文は `Proposition 1.4, (iv)` から取る。

  原文 (GenEll p.17):
  > Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

  ★`ht^Falt` が有界な点は有限個しかない。 -/
  northcott : ∀ (C : ℝ) (d : ℕ), 0 < d → {x ∈ degLe d | faltingsHeight x ≤ C}.Finite
  /-- 原文 `E′ ≙ E/H_F`——l-cyclic 部分群スキームによる商。

  原文 (GenEll p.17):
  > Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

  ★`Lemma 3.5` の証明はこの商を取って `Proposition 3.4` を適用する。 -/
  quotLCyclic : Curve → ℕ → Curve
  /-- ★★★★★★**`deg_∞(E′) = l·deg_∞(E)`**——`Lemma 3.2` の**大域版**。

  原文 (GenEll p.17):
  > Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

  ★局所版は `Interface/GenEll/TateLocal.lean` の上の `Lemma 3.2, (ii)` であり、
  そこから素点にわたって足し上げる段は原文に書かれていない。
  ★★`l` が局所高さと素であることがここで効く
  ——`Lemma 3.2, (i)` により `H` は各素点で `μ_l` になる。 -/
  degInf_quotLCyclic : ∀ (E : Curve) (l : ℕ), Nat.Prime l → HasLCyclic E l →
    PrimeToLocalHeights E l →
    degInf (cls (quotLCyclic E l)) = (l : ℝ) * degInf (cls E)
  /-- ★★★★★**`ht^Falt(E′) ≤ ht^Falt(E) + 2log(l) + C₀`**——l-同種による変化。

  原文 (GenEll p.17):
  > Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

  ★原文の典拠は **[FC] Chapter I, Proposition 2.7**(次数 `l` の被覆射の延長)と、
  「(1,1)-形式を `E_v` 上で積分するのと `(E_H)_v` 上で積分するのとで `l` 倍だけ違う」
  という複素解析の段である(原文は 1 文で済ませている)。
  ★★**`C₀` は `E`･`l` に依らない**——だから `∃` を外側に置く。 -/
  faltingsHeight_quotLCyclic : ∃ C₀ : ℝ, ∀ (E : Curve) (l : ℕ), Nat.Prime l → HasLCyclic E l →
    faltingsHeight (cls (quotLCyclic E l))
      ≤ faltingsHeight (cls E) + 2 * Real.log l + C₀
  /-- ★`d = [L:ℚ]` は正。 -/
  degOfDefinition_pos : ∀ E : Curve, 0 < degOfDefinition E
  /-- ★★★★★**`l` が十分大きければ局所高さと素になる**。

  原文 (GenEll p.18):
  > Lemma 3.7. (Finite Exceptional Sets) Let

  ★原文は『if `v` is any local height of `E_L`, then `d·deg_∞([E_L]) ≥ v·log(2)`』を
  証明なしで使う。★★`l` が素数ですべての局所高さより大きければ、
  `l` はそれらを割らない——この 2 段をまとめて欄に出している。 -/
  primeToLocalHeights_of_lt : ∀ (E : Curve) (l : ℕ), Nat.Prime l → SemiStable E →
    (degOfDefinition E : ℝ) * degInf (cls E) < (l : ℝ) * Real.log 2 →
    PrimeToLocalHeights E l
  /-- ★★★★★★**l-cyclic かつ局所高さと素な類の例外集合**。

  原文 (GenEll p.18):
  > Lemma 3.7. (Finite Exceptional Sets) Let

  ★これが Galois-finite であることは **`Lemma 3.5`**(高さの不等式)と
  **`Lemma 3.6`**(初等的な評価)から `ht^Falt` が有界になり、
  **`Proposition 1.4, (iv)`**(Northcott)と **`Example 1.3, (i)`** で出る。 -/
  lcyclicExc : Set EllClass
  galoisFinite_lcyclicExc : GaloisFinite lcyclicExc
  mem_lcyclicExc : ∀ (E : Curve) (l : ℕ), Nat.Prime l → SemiStable E →
    HasLCyclic E l → PrimeToLocalHeights E l → cls E ∈ lcyclicExc
  /-- ★★★★★**compactly bounded の中で乗法還元を持たない類の例外集合**。

  原文 (GenEll p.18):
  > Lemma 3.7. (Finite Exceptional Sets) Let

  ★compactly bounded な集合(`Example 1.3, (ii)`)の中で、
  乗法還元を持たないものは潜在的に良還元であり、高さが有界になる。 -/
  noMultRedExc : Set EllClass → Set EllClass
  galoisFinite_noMultRedExc : ∀ KV : Set EllClass, CompactlyBounded KV →
    GaloisFinite (noMultRedExc KV)
  mem_noMultRedExc : ∀ (KV : Set EllClass) (E : Curve), cls E ∈ KV → ¬ HasMultRed E →
    cls E ∈ noMultRedExc KV
  /-- ★Galois-finite は有限合併で閉じる(`Example 1.3, (i)`)。 -/
  galoisFinite_union : ∀ S T : Set EllClass, GaloisFinite S → GaloisFinite T →
    GaloisFinite (S ∪ T)

/-- ★Track B は何を作らねばならないか。 -/
def EllModuliData.waiting : WaitingFor :=
  { what := "楕円曲線のモジュライスタック M̄_ell と無限遠因子、Faltings 高さ ht^Falt、deg_∞ / ht_∞、log-diff_{M̄_ell}、および半安定還元・l-cyclic 部分群スキーム・最小定義体"
    trackB := "Found/GenEll — ★ht^Falt は Arakelov 理論(§1)を、SemiStable / HasLCyclic は Tate 曲線(Interface/GenEll/TateLocal.lean)を要求する。★log-diff は Definition 1.5, (iii) すなわち ADiv(F) と degNormalized の上に立ち、その 2 つは Found/GenEll/ArithDiv.lean に実装済みである" }

/-! ## ★出典の紐付け(`.src`) -/

def EllModuliData.src : Source :=
  { paper := "GenEll", pdfPage := 17, item := "Proposition 3.4",
    sectionId := "genell-prop-3-4" }

end ABC3.Interface.GenEll
