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
  /-- 原文「`E_L` admits an **l-cyclic** subgroup scheme `H_L ⊆ E_L`」。

  ## ★★★★★★2026-08-29 の測定 —— **標数 0 では「Galois 安定な直線」と同じ**

  ★原文の `H_F ⊆ E_F` は**生成ファイバー**（数体 `F` の上）の部分群スキームである。
  ★★標数 0 では有限群スキームはすべてエタール（Cartier）であり、
  エタール有限群スキームは**有限 Galois 加群**と圏同値である。
  ★★★したがって

      「`l`-巡回部分群スキーム `H_F ⊆ E_F`」 ⟺ 「`E[l]` の中の `Gal(ℚ̄/F)`-安定な直線」

  である。**有限平坦群スキームの一般論は要らない**
  （それが要るのは `Lemma 3.2, (ii)` の `E/μ_l` を `𝒪_K` 上で作る段だけである）。

  ★★★★これで `Theorem 3.8` の側では、`¬ HasLCyclic` は
  「安定な直線が無い」と読め、`Found/GenEll/Thm38Bridge.lean` の
  `exists_nonUpper_of_no_stable_line` がそのまま効く。 -/
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
  ——`Lemma 3.2, (i)` により `H` は各素点で `μ_l` になる。

  ★★★**2026-09-02（第 1340）の締め直し**——かつては**等式**であったが、
  消費側（`Section3.lean` の `lemma_3_5`）が使うのは
  `l·deg∞(E) ≤ deg∞(E′)` の**向きだけ**である。
  ☆こう弱めると、良い素点では `0 ≤ minDeltaExp` が自明になり、
  この欄は**悪い素点の Tate の関係だけ**（在庫、第 1141）に落ちる。 -/
  degInf_quotLCyclic : ∀ (E : Curve) (l : ℕ), Nat.Prime l → HasLCyclic E l →
    PrimeToLocalHeights E l →
    (l : ℝ) * degInf (cls E) ≤ degInf (cls (quotLCyclic E l))
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
  **`Proposition 1.4, (iv)`**(Northcott)と **`Example 1.3, (i)`** で出る。

  ## ★★★★★★★★2026-08-31 の訂正(第 755)——`l` の下界が要る

  ☆以前ここは

      lcyclicExc : Set EllClass
      mem_lcyclicExc : ∀ E l, Prime l → SemiStable E →
        HasLCyclic E l → PrimeToLocalHeights E l → cls E ∈ lcyclicExc

  であった。★★**これは強すぎて witness が作れない**
  (`Check/GenEll/LcyclicExcTooStrong.lean` の測定):
  `Lemma 3.5` が与えるのは `(l/14)·deg∞ ≤ ht^Falt + 2·log(l) + C′` であり、
  `ht^Falt` の**上界**を出すには `l` が `ht^Falt` に比べて大きいこと
  ——原文の**条件 (a)**——が要る。`l = 2` を取れば対象は無限個になる。

  ★★★`Skeleton/GenEll/Section3.lean` の `lemma_3_7` は `mem_lcyclicExc` を呼ぶ時点で
  `condA ∨ condB` を持っているので、**それを仮説として渡す**形に直した。
  ☆`lemma_3_7`・`theorem_3_8`・`Corollary 4.3/4.4` の**statement は変わらない**。 -/
  lcyclicExc : ℝ → ℝ → Set EllClass → Set EllClass
  galoisFinite_lcyclicExc : ∀ (C eps : ℝ) (KV : Set EllClass), CompactlyBounded KV →
    GaloisFinite (lcyclicExc C eps KV)
  mem_lcyclicExc : ∀ (C eps : ℝ) (KV : Set EllClass) (E : Curve) (l : ℕ),
    Nat.Prime l → SemiStable E → HasLCyclic E l → PrimeToLocalHeights E l →
    (((100 * (degOfDefinition E : ℝ)
          * (faltingsHeight (cls E) + C * (degOfDefinition E : ℝ) ^ eps) ≤ (l : ℝ))
        ∧ HasMultRed E)
      ∨ cls E ∈ KV) →
    cls E ∈ lcyclicExc C eps KV
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
  /-- ★★★★★★**3･5 捧れを有理化する拡大 `L′`**(`Theorem 3.8` の証明の骨格)。

  原文 (GenEll p.19):
  > Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

  ★原文 p.20:『there exists a Galois extension `L′` of `L` of degree that divides
  `d₀ = (3²−1)(3²−3)(5²−1)(5²−5) = 23040`、so as to render the 3- and 5-torsion
  points of `E_L` rational over `L′` … we may assume that `E_{L′}` has **semi-stable
  reduction** at all of the finite primes』。 -/
  torsionExt : Curve → Curve
  /-- ★`ℚ̄` 上の同型類は変わらない——だから `ht^Falt` も `deg_∞` もそのまま。 -/
  cls_torsionExt : ∀ E : Curve, cls (torsionExt E) = cls E
  /-- ★★次数は高々 `23040` 倍——これが `Theorem 3.8` の係数 `23040` の出所である。 -/
  degOfDefinition_torsionExt : ∀ E : Curve,
    degOfDefinition (torsionExt E) ≤ 23040 * degOfDefinition E
  /-- ★★★**`L′` の上では半安定還元になる**——原文の『we may assume』の中身。 -/
  semiStable_torsionExt : ∀ E : Curve, SemiStable (torsionExt E)
  /-- ★★★**潜在的乗法還元は `L′` 上で乗法還元になる**。 -/
  hasMultRed_torsionExt : ∀ E : Curve, HasPotMultRed E → HasMultRed (torsionExt E)
  /-- ★★★★★**`30` と互いに素なら局所高さと素であることは `L′` へ移る**。

  ★原文 p.20 の括弧:『passing to such a Galois extension of `L` only affects the prime
  decomposition of the local heights via the primes that divide `d₀`, of which there are
  only finitely many, namely, **2, 3, and 5**』。
  ★★これが `Theorem 3.8` の条件 (b) の『`30 = 2·3·5` と素』の出所である。 -/
  primeToLocalHeights_torsionExt : ∀ (E : Curve) (l : ℕ), PrimeToLocalHeights E l →
    Nat.Coprime l 30 → PrimeToLocalHeights (torsionExt E) l
  /-- ★★★★★★**原文 p.20 の最終段**。

  原文 (GenEll p.19):
  > Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

  ★乗法還元の素点で局所高さが `l` で割れなければ、局所理論(`Lemma 3.2` の直前)により
  mod `l` 像は `α = (1 1 / 0 1)` を含む。
  ★★`l`-巡回部分群スキームを持たなければ、mod `l` 像は**非上三角**行列を含む。
  ★★★その 2 つから **`Lemma 3.1, (iv)`** で `GL₂(ℤ_l)` の像は `SL₂(ℤ_l)` を含む。
  ★★★★`Lemma 3.1` は `Found/GenEll/Lemma31.lean`･`Sl2Padic.lean` に**4 条すべて実装済み**である。

  ## ★★★★★★★★★★ 2026-08-29 の再実測——**前の記述は古かった**

  ★以前ここには「**Galois 表現そのものが未構築**なので」と書いてあったが、
  ★★**それはもう当てはまらない**——`Found/GaloisRep/GalRep.lean` に

      `galRep : Gal(L/K) →* GL₂(ℤ_l)`（`exists_galRep` もある）

  が構成されている。★`Found/GaloisRep/` は 325 ファイルあり、
  Tate 加群・Weil 対・行列式＝円分指標まで入っている。

  ★★★**したがってこのフィールドが本当に受けているのは次の 2 つだけ**である:

  1. 局所理論の行列表示——`0 → F_l(1) → M_l(E) → F_l → 0` に合わせた基底で
     `l ∤ v(q)` なら mod `l` 像が `α = (1 1 / 0 1)` を含むこと。
     ★Tate 一意化と `Lemma 3.2, (i)` は**閉じている**（`Lemma32Tate.lean`）。
  2. `l`-巡回部分群スキームを持たないこと ⟹ mod `l` 像が**非上三角**行列を含むこと
     （`l`-巡回 ⟷ Galois 安定な直線の対応）。

  ★★★★★★★**`Theorem 3.8` に Serre の開像定理は要らない**——群論の核は
  `Lemma 3.1, (iv)` だけであり、それは済んでいる。
  ★★★★★結論が `E`(`L` 上)についてなのは、
  `Im(Gal_{L′}) ⊆ Im(Gal_L)` だからである。

  ## ★★★★★★★★2026-08-31 の訂正(第 776)——`5 ≤ l` が要る

  ☆以前ここに `5 ≤ l` は無かった。★★**群論の核 `Lemma 3.1, (iv)` が `5 ≤ l` を
  要求する**(`Found/GenEll/Sl2Padic.lean` の `lemma_3_1_iv`)ので、
  `l ∈ {2, 3}` では witness が作れない。
  ★`Check/GenEll/ImageSL2NeedsL5.lean` の測定を参照。

  ★★★`Theorem 3.8` の側では両方の条件から `5 ≤ l` が出るので、statement は変わらない:
  条件 (a) では `l ≥ 23040·100 ≥ 5`、条件 (b) では `l` は `30` と素な素数だから `l ≥ 7`。 -/
  imageContainsSL2_of_torsionExt : ∀ (E : Curve) (l : ℕ), Nat.Prime l → 5 ≤ l →
    HasMultRed (torsionExt E) → PrimeToLocalHeights (torsionExt E) l →
    ¬ HasLCyclic (torsionExt E) l → ImageContainsSL2 E l
  /-- ★★★★★**`SL₂` を含むことと全射性は、`l` が `L` で不分岐なら同じこと**である。

  原文 (GenEll p.22):
  > Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

  ★原文 p.22 の冒頭:『if `l` is any prime number that is **unramified in `L`**, then the
  image of the Galois representation … contains `SL₂(ℤ_l)` **if and only if** the Galois
  representation … is **surjective**』。
  ★★理由も原文が括弧で書いている——`ℚ(ζ_{l^∞})/ℚ` は `l` で**完全分岐**するので
  `L/ℚ` と**線型無関連**であり、円分円指標が全射になる。
  ★★★使うのは片向きだけなので、そちらだけを欄に出す。 -/
  imageSurjective_of_containsSL2 : ∀ (E : Curve) (l : ℕ), Nat.Prime l →
    PrimeToRamification E l → ImageContainsSL2 E l → ImageSurjective E l
  /-- ★空集合は compactly bounded である(`Example 1.3, (ii)`)。

  ★`Corollary 4.3` は `K_V` を持たないが、`Theorem 3.8` の条件 (a) を使うには
  何か 1 つ compactly bounded な集合を渡す必要がある——条件 (a) は `K_V` を使わない。 -/
  compactlyBounded_empty : CompactlyBounded (∅ : Set EllClass)
  -- ★★★★§4 —— 原文 p.22-23 が Corollary 4.3 の証明で引くもの
  /-- ★★★★**乗法還元の素点の下にある `ℚ` の素数と、そこでの局所高さ**。

  原文 (GenEll p.22):
  > Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

  ★原文 p.23:『the “`h`” of Lemma 4.2 corresponds to **`23040d · deg∞([E_L])`**
  [cf. the meaning of “`d₀ = 23040`” in the proof of Theorem 3.8]』。
  ★★すなわち `∑ h_j·log(p_j) = 23040·d·deg∞` であり、これが `Lemma 4.2` への入力になる。 -/
  multCard : Curve → ℕ
  multCard_pos : ∀ E : Curve, 0 < multCard E
  multPrime : ∀ E : Curve, Fin (multCard E) → ℕ
  multPrime_prime : ∀ (E : Curve) (j : Fin (multCard E)), Nat.Prime (multPrime E j)
  localHt : ∀ E : Curve, Fin (multCard E) → ℕ
  localHt_pos : ∀ (E : Curve) (j : Fin (multCard E)), 0 < localHt E j
  /-- ★★★★★**`∑ h_j·log(p_j) = 23040·d·deg∞([E_L])`**(原文 p.23)。 -/
  sum_localHt_eq : ∀ E : Curve,
    (∑ j : Fin (multCard E), (localHt E j : ℝ) * Real.log (multPrime E j))
      = 23040 * (degOfDefinition E : ℝ) * degInf (cls E)
  /-- ★★★**`S∘` のうち `S` 以外の部分**。

  原文 p.22:『write `S∘` for the union of `S`, the primes of `ℚ` that lie under primes of
  **potentially multiplicative reduction** of `E_L`, and the primes that appear in the
  **prime decomposition of the local heights** of `E_L`』。 -/
  badPrimes : Curve → Finset ℕ
  badPrimes_prime : ∀ (E : Curve), ∀ p ∈ badPrimes E, Nat.Prime p
  /-- ★★★★**`Lemma 4.2` が押さえるのはこの和である**。

  ★左辺は `x_{S∘} − x_S`、右辺の第 1 項が乗法還元の素数、第 2 項が
  局所高さ `h_j` の素因数分解に現れる素数である。
  ★★この形にしておけば、**実装済みの `Lemma 4.2`** がそのまま使える。 -/
  sum_log_badPrimes_le : ∀ E : Curve, (∑ p ∈ badPrimes E, Real.log p)
      ≤ (∑ j : Fin (multCard E), Real.log (multPrime E j))
        + (∑ j : Fin (multCard E), Real.log ((localHt E j : ℝ) + 1))
  /-- ★★`badPrimes` の外の素数は条件 (a) の前半を満たす。 -/
  primeTo_badPrimes : ∀ (E : Curve) (l : ℕ), Nat.Prime l → l ∉ badPrimes E →
    PrimeToMultPrimes E l ∧ PrimeToLocalHeights E l
  /-- ★★★**`S•` のうち `S` 以外の部分**。

  原文 p.22:『write `S•` for the union of `S∘`, the primes of `ℚ` that **ramify in `L`**,
  and the primes that divide the **ramification indices** of primes of `ℚ` in `L`』。 -/
  ramPrimes : Curve → Finset ℕ
  ramPrimes_prime : ∀ (E : Curve), ∀ p ∈ ramPrimes E, Nat.Prime p
  badPrimes_subset_ramPrimes : ∀ E : Curve, badPrimes E ⊆ ramPrimes E
  primeTo_ramPrimes : ∀ (E : Curve) (l : ℕ), Nat.Prime l → l ∉ ramPrimes E →
    PrimeToRamification E l
  /-- ★★★★**`x_{S•} ≤ x_{S∘} + 3d·log-diff`**。

  原文 p.23:『since [as is easily verified, by considering the **trace** of an extension of
  number fields] the primes appearing in the arithmetic divisor that gives rise to
  “`log-diff_Mell`” [cf. Definition 1.5, (iii)] appear with **multiplicity ≥ one less than
  the ramification indices** of `L/ℚ`』。
  ★これが `Corollary 4.3, (c)` の `6d·log-diff` の出所である(2 倍される)。 -/
  sum_log_ramPrimes_le : ∀ E : Curve, (∑ p ∈ ramPrimes E, Real.log p)
      ≤ (∑ p ∈ badPrimes E, Real.log p)
        + 3 * (degOfDefinition E : ℝ) * logDiffMell (cls E)

/-- ★Track B は何を作らねばならないか。 -/
def EllModuliData.waiting : WaitingFor :=
  { what := "楕円曲線のモジュライスタック M̄_ell と無限遠因子、Faltings 高さ ht^Falt、deg_∞ / ht_∞、log-diff_{M̄_ell}、および半安定還元・l-cyclic 部分群スキーム・最小定義体"
    trackB := "Found/GenEll — ★ht^Falt は Arakelov 理論(§1)を、SemiStable / HasLCyclic は Tate 曲線(Interface/GenEll/TateLocal.lean)を要求する。★log-diff は Definition 1.5, (iii) すなわち ADiv(F) と degNormalized の上に立ち、その 2 つは Found/GenEll/ArithDiv.lean に実装済みである" }

/-! ## ★出典の紐付け(`.src`) -/

def EllModuliData.src : Source :=
  { paper := "GenEll", pdfPage := 17, item := "Proposition 3.4",
    sectionId := "genell-prop-3-4" }

end ABC3.Interface.GenEll
