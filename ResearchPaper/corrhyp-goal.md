# [CorrHyp] トラックのゴール(2026-09-04 設定、ABC3b)

対象: S. Mochizuki, *Correspondences on Hyperbolic Curves* [CorrHyp]
(物理 18 ページ、著者 1 名、`papers.json` に短縮タグ登記済み、`0_Source` に PDF/txt あり)。

**★2026-09-04、p.3-13・p.15・p.17 を 200dpi 目視確認済み**(`papers.json`
`notationRisk: "medium"`)。実害は overline(X̄ → X)と script 文字(𝒳・𝒴 → X・Y)の
脱落の 2 点(GenEll のような行列の順序入れ替えは無い、詳細は `papers.json` を見よ)。

再現: `node tools/paper-items.mjs CorrHyp`(ABC3c が LocProP のために新設した汎用ツール。
行頭 `Kind N.M.`(ピリオド+行末 or 開き括弧)の宣言規則で番号付き項目を §ごとに数える)。
24 件・重複なしを確認済み(手作業での逐次読解と突き合わせて一致)。

---

## 0. ゴール(現在地)

> **CorrHyp §1 5/5, §2 6/6, §3 3/3, §4 2/2, §5 7/7, §6 1/1 —— 合計 24/24(Skeleton)**

★★**2026-09-04、S0(Skeleton = statement を型で固定する段)を完了した。**
`Interface/CorrHyp/HyperbolicCurve.lean`(`HyperbolicCurveData`・`StackType` を posit)
と `Skeleton/CorrHyp/Section1.lean`〜`Section6.lean` に 24 項目すべてを置き、
`lake build` 0 エラー・`tools/check.mjs` で G1(出典・逐語照合)を全項目パスさせた。
残る `sorry` は 19 件(Definition 5 件は sorry 無しで構成済み)。G9(非空虚性の対照)は
14 件が未着手のまま——プロジェクト全体でトラッキング対象の既知 debt であり、
`check.mjs` 自身が「新規は落とす(ブロックしない)」と明記している。
次の段(S1: Track B で `Found/CorrHyp/` を積む)はまだ手つかず。

## 1. §ごとの内訳(節タイトルつき)

| § | 節タイトル | 件数 | 内訳(宣言順) | 物理ページ |
|---|---|---:|---|---|
| 1 | Basic Definitions | 5 | Definition 1.1 / 1.2 / 1.3 / 1.4 / 1.5 | 3–4 |
| 2 | Review of Results of Margulis and Takeuchi | 6 | Definition 2.1 / 2.2 / 2.3, Proposition 2.4, Theorem 2.5 / 2.6 | 5–7 |
| 3 | The Non-arithmetic Case | 3 | Definition 3.1, Proposition 3.2, Theorem 3.3 | 8 |
| 4 | The Main Theorem | 2 | Lemma 4.1, Theorem 4.2 | 9–10 |
| 5 | Isogenies of General Curves | 7 | Lemma 5.1, Definition 5.2, Theorem 5.3, Lemma 5.4 / 5.5 / 5.6, Theorem 5.7 | 11–15 |
| 6 | Interpretation of a Theorem of Royden | 1 | Theorem 6.1 | 17 |

★`Corollary` は本論文に 0 件。`Remark` は本論文では**すべて無番号**(`"Remark."` のみ、
`Remark N.M` 形式が無い)ので LocProP と同様に項目数に入らない。

★§0(Introduction)は分母 24 から**除外した**(LocProP の §0 とは扱いが異なる——
LocProP の §0 は新規の Definition/Lemma を持つが、本論文の §0 は `Theorem A`/`B`/`C` の
3 件しかなく、そのいずれも本文中で明示的に後続定理の再掲と言明されている。二重計上を
避けるための除外であり、GenEll の §1 起点と同じ扱い):

| §0 の呼称 | 本文中の再掲先の明言(逐語) | ページ |
|---|---|---|
| Theorem A | 「cf. Lemma 4.1 and Theorem 4.2 in the text」 | p.1 |
| Theorem B | 「Theorem B follows from Theorem 5.3 in the text」 | p.2 |
| Theorem C | 「(given as Theorem 6.1 in the text)」 | p.2 |

## 2. 主定理との対応(★本文の言明そのもの、LocProP と違い推定ではない)

| 呼称 | 内容(導入部の記述) | 対応する節 |
|---|---|---|
| Theorem A | 双曲曲線に isogenous な曲線の有限性(第一主定理) | `Lemma 4.1` + `Theorem 4.2`(§4) |
| Theorem B | 一般の曲線は自身の hyperbolic core に一致——非自明な相関を持たない | `Theorem 5.3`(§5) |
| Theorem C | `M_{g,r}` は非自明な自己同型・相関を持たない(Royden の定理の帰結) | `Theorem 6.1`(§6) |

§2(Margulis・Takeuchi の結果の review、`Theorem 2.5`/`2.6`)は他論文 `[Marg]`/`[Take]` の
引用結果であり、本プロジェクトでは証明せず posit する対象になる可能性が高い(要検討——
`Interface` 行き。`genell-track-b` で `Theorem 3.8` 等の外部結果を扱った前例を参照)。

## 3. Track B(Found/CorrHyp/)——2026-09-04 に着手、実測に基づく規模感

★★**Found(sorry 無しの証明)まで24件全部を1セッションで終える見込みは無い**
——GenEll(数か月・現在も§3-§4に3本のsorryが残る)と同型の規模になることを、
実際に mathlib を調べたうえで確認した。「壁」ではなく道として、葉から実装する。

### 使える資産(mathlib、2026-09-04 実測)

| 材料 | 場所 |
|---|---|
| `ℍ`(上半平面)・`SL(2,ℝ)` の Möbius 作用 | `Analysis.Complex.UpperHalfPlane.*` |
| `Subgroup 𝒢 ≤ SL(2,ℝ)`(離散)→ `ℍ` に固有不連続作用 | `UpperHalfPlane.instProperlyDiscontinuousSL2RSubgroup` |
| Schreier の補題(有限指数部分群は有限生成) | `Subgroup.fg_of_index_ne_zero` |
| 被覆空間・モノドロミー | `Topology.Covering.*` |
| `PSL(n,K)` の構成 | `Matrix.ProjectiveSpecialLinearGroup` |

### 無い(§2・§6 の律速)

Margulis/Shimura-arithmeticity(代数群・非可換 Galois コホモロジー・Brauer 群)、
四元数環の分類理論、モジュライスタック `M_{g,r}`、Teichmüller 空間、
Royden の定理——いずれも GenEll のモジュラー多項式ギャップと同型の
「未構築の数学」で、§2(`Definition 2.2`/`2.3`・`Proposition 2.4`・`Theorem 2.5`-`2.6`)
と §6(`Theorem 6.1`)はこれ単体で人年規模になりうる。

### 実装した分(`Found/CorrHyp/FuchsianGroup.lean`、すべて sorry 無し・標準3公理のみ)

★2026-09-04、`Space := FuchsianGroup`(`SL(2,ℝ)` の離散部分群、原文の
`PSL₂(ℝ)⁰` からの逸脱を記録)というモデルで §1・§3 の骨格を実装した:

| 節点 | 実装した宣言 | 状態 |
|---|---|---|
| §1 `Definition 1.1`-`1.4` の土台 | `FEt := IsFiniteIndexIn`・`comp`・`pullback`(`inter`)・`pbFst`・`pbSnd`・`idFEt` | ✅ 圏構造が揃った |
| §1 `Definition 1.5` 直後(isogeny は同値関係) | `isIsogenous_refl`/`_symm`/`_trans` | ✅ **完成** |
| §2 `Definition 2.1` 直前(`Γ ⊆ Comm(Γ)`) | `self_le_commensurator`(mathlib の `Subgroup.Commensurable.commensurator` を発見・流用) | ✅ **完成** |
| §3 `Proposition 3.2`(`Γ_Z ⊆ Comm(Γ_X)`) | `prop_3_2`(`Ncore` = `Z.Γ` の中での `C.Γ` の normal core、を経由) | ✅ **完成**——numbered item の本体として初 |
| §3 `Theorem 3.3` | `fg_of_finiteIndexIn`(Schreier の補題、`Theorem 3.3` が要る片方向) | 部分——下記参照 |

### `Theorem 3.3` に残っていた 2 つの壁

1. **群論側 —— ★2026-09-04 解消**: 原文は「`Comm(Γ_X)` は有限生成部分群 `Γ_X`
   を有限指数に持つので、それ自身も有限生成」と1行で済ませているが、これは
   `fg_of_finiteIndexIn`(有限指数**部分群**が有限生成、Schreier)の**逆向き**
   ——「有限生成な有限指数部分群を持つ群は有限生成」——であり、`exact?` で
   mathlib に見当たらなかった。`fg_of_fg_finiteIndex`(Reidemeister–Schreier の
   最も単純な形、`H` の生成元 + `G/H` の代表系が `G` を生成する)として
   自分で証明した(sorry 無し)。
2. **幾何側(まだ未着手、規模を実測した)**: 「次数が `g',r'` と `e_Y` だけで
   抑えられる」という原文の主張は、`FuchsianGroup` を実際の種数・穴の数
   (§1 冒頭の `type`)に結び付ける Riemann–Hurwitz 型の計算を要求する。
   ★★**mathlib に Gauss–Bonnet・双曲多様体の面積・Riemann 面の種数の理論は
   0 件**(2026-09-04、`grep` で `GaussBonnet`/`RiemannSurface`/幾何的な
   `genus` を実測——`EulerCharacteristic` はホモロジー代数の鎖複体のものだけ)。
   これは GenEll のモジュラー多項式ギャップと同型の「未構築の数学」であり、
   `§3 単体でもTheorem 3.3・Theorem 4.2・Theorem 5.3・Lemma 5.4`-`5.7` 全体の
   共通の律速——**Riemann 面の双曲幾何(面積・種数・Gauss–Bonnet)を
   mathlib レベルで新たに構築する必要がある人年規模の欠落**として記録する。
   `Skeleton/CorrHyp/Section5.lean` の `StackType`/`e_Y` はこの欠落を
   posit で受けている(Skeleton 段階で正しい選択だった)。

### 次の一手(依存の少なさ順)

★★純粋な群論(`FuchsianGroup` の上だけで閉じる節点)は §1・§3 の
`Proposition 3.2` まで**掘り尽くした**——残る節点はどれも「未構築の数学」
(幾何・スキーム論・代数群論)を新たに要求する側に移った:

1. §5 の `StackType`/`e_Y`(Euler 標数)を `FuchsianGroup` の商 `ℍ/Γ` の
   実際の双曲幾何(面積・種数)と結び付ける——mathlib に Gauss–Bonnet が
   無いので、これ自体が新規のライブラリ構築になる。できて初めて
   `Theorem 3.3`/`Theorem 4.2`/`Theorem 5.3`/`Lemma 5.4`-`5.7` に届く。
2. §4(`Lemma 4.1`・`Theorem 4.2`)は係数拡大 `Ext`(k → K)が要り、
   スキーム論(有限型・降下理論)を要求する——①とは別種の未構築の数学。
3. §2 の `MargulisArithmetic`/`ShimuraArithmetic` 本体と §6 の Royden の定理は
   代数群論・Teichmüller 空間という、①②とも異なる第三の未構築の数学を要する。
4. G9(非空虚性の対照)は Track B が本物のデータを供給して初めて
   意味のある形で作れる(`.waiting` に記録済み)。

★どの道も本質的に異なる分野の「未構築の数学」に当たる——GenEll が
モジュラー多項式ギャップ1つで止まったのと違い、本論文は最初から
3つの独立した大きな塊(双曲幾何・スキーム論・代数群論)を要求する
構造になっている(2026-09-04 実測、上記①②③)。

★★**3つの塊のうち②(スキーム論)だけ、mathlib に相応の蓄積がある**
(2026-09-04 追加実測)。`AlgebraicGeometry.Scheme.Etale`・
`LocallyOfFiniteType` に加えて `AlgebraicGeometry/SpreadingOut.lean` に
「有限型なら有限生成部分環まで下りて一意に持ち上がる」という
spreading-out の補題群が存在する。①(双曲幾何、Gauss–Bonnet 0件)・
③(代数群論、非可換 Galois コホモロジー)より見込みがある。

★★★★★★★**2026-09-04、セッション末に①を再実測・確定**
(`genus`・`GaussBonnet`・`RiemannHurwitz` で `grep`——3語とも**文字通り
0件**、無関係な `HurwitzZeta`(数論的ゼータ関数、別分野)のみヒット)。
これで「双曲幾何(種数・Euler標数・Gauss–Bonnet)が mathlib に無い」は
複数回の実測すべてで確認された確定事実——`Theorem 3.3`・§5 の大半・
間接的に `Theorem 4.2`/`Theorem 6.1` の律速。次にここへ戻るセッションは
再確認不要、この人年規模の欠落を埋める作業(mathlib への新規貢献に近い
規模)から直接始めてよい。

★★★**さらに一段実測(2026-09-04)——`Lemma 4.1` に直撃する道具を発見**:
`AlgebraicGeometry/AffineTransitionLimit.lean` は「スキームが affine な
遷移射を持つ余濾過的な図式(`D : I ⥤ Scheme`)の極限である」ときの
spreading-out 定理群を持つ。これは `Spec K`(`K` を有限生成 `k`-部分環
`R_i ⊆ K` の余極限とみなしたときの極限スキーム)の設定そのものであり、
`Lemma 4.1` が要る2つの向きにそれぞれ対応する補題が**両方存在する**:

| `Lemma 4.1` が要ること | 対応する mathlib 補題 |
|---|---|
| `K` 上で定義された対象は、ある有限段 `R_i` まで下りて存在する(spreading out) | `Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation` |
| 極限上で一致する2つの射は、ある有限段まで下りればすでに一致する(rigidity/uniqueness) | `Scheme.exists_hom_comp_eq_comp_of_locallyOfFiniteType` |

残る作業は「`X_K`(`K`上の双曲曲線)」「`Z_K` との correspondence」を
この抽象的な極限の枠組みに具体的に載せる段——`I` を `K` の有限生成
`k`-部分環の余濾過圏に取り、`D` をそれぞれの `Spec R_i` 上の底変換の
図式にする、という設定作業がまだ残っている(未着手、規模は未測定だが
道具はすでに mathlib にある)。

★★★★**2026-09-04続報(第23-24件)、上の「設定作業」の骨組みを完成**。
`Found/CorrHyp/FieldLimit.lean`:
- `I := (FgSubalgebra k K)ᵒᵖ`(`K` の有限生成 `k`-部分環全体、包含で有向)
  を実際に構成し、`D := toSchemeDiagram k K := (toRingCat k K).op ⋙ Scheme.Spec`
  という図式を作った。
- **`Spec K` が `D` の極限であること**(`isLimit_specKCone`)を sorry 無しで
  証明した——環側で `K` がその有限生成部分環の余極限であること
  (`isColimitToRingCatCocone`、`Subalgebra.iSupLift` の `RingHom` 版が
  mathlib に無かったので手で構成)を示し、`IsColimit.op` + `Scheme.Spec` が
  極限を保つこと(`Γ ⊣ Spec` の右随伴、`Adjunction.rightAdjoint_preservesLimits`)
  で `Scheme` 側へ運んだ。
- 上の表の2定理が要求する側条件(`IsAffineHom(D.map f)`・
  `CompactSpace(D.obj i)`・`QuasiSeparatedSpace(D.obj i)`)を、
  `D.obj`/`D.map` が常に `Spec(…)` の形であることから**無条件の instance**
  として確立した(mathlib の一般事実だけで済んだ)。

**結果**: `c := specKCone k K`・`hc := isLimit_specKCone k K` として、
`Scheme.exists_hom_comp_eq_comp_of_locallyOfFiniteType`/
`Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation` を**側条件抜きで
直接呼べる状態**になった——`Lemma 4.1` のスキーム論的な骨組みは完成。

**まだ残る(誠実な現状)**: `lemma_4_1`(Skeleton)自体は `HyperbolicCurveData`
の抽象フィールド(`Space`・`Ext`・`IsGenericallyScheme`・`ModuliStack`)で
書かれており、これらを「`Spec K` 上の実際のスキーム」に結び付ける
**橋渡しフィールド**が `Interface` にまだ無い。橋渡しを足しても、
`.needs` に記録済みの `LocallyOfFinitePresentation f`(`X_K → Spec k`)や
接空間の計算(標数0での単射性)は依然として要る——モジュライスタック
`M_{g,r}` の実現(§5、mathlib に不在)と分かちがたく結びついている。
現時点で numbered item として Found 完成しているのは
`Proposition 3.2`(§3)のみで、24/24 にはまだ遠い。

★★★★★**2026-09-04続報、`Comm : Fuchsian → Fuchsian` の型付けが
Margulis の二分法そのものを埋め込んでいることを原文照合で確認**(重要な発見)。
`corrhyp_prop_3_2`(`Found/CorrHyp/FuchsianGroup.lean`)に `.src` を付けて
Skeleton の `prop_3_2` と直結させようとして、次に気づいた:

- 原文 p.5(`Comm(Γ)` 導入直後): 「X is not arithmetic. Thus, we have
  Γ ⊆ Comm(Γ) ⊆ PSL2(R)⁰, and **Γ is of finite index in Comm(Γ)**.」
  ——非 arithmetic ⟹ Comm(Γ) 自身が離散(Fuchsian)、という**Margulis の
  二分法の主張そのもの**が §2 の地の文で明言されている。
- `Interface/CorrHyp/HyperbolicCurve.lean` の `Comm : Fuchsian → Fuchsian`
  という型(`Fuchsian` = 離散部分群)は、この事実を**型シグネチャの中に
  既に組み込んでいる**——`Definition 2.1`(`InfinitelyManyCorr`)・
  `Proposition 3.2` の**両方**が `D.Comm Γ : D.Fuchsian`(離散である前提)を
  直接使っているため、これは Skeleton 設計時の見落としではなく、
  原文の記述を素直に型にした結果。
- ★しかし `Comm` が**全域関数**である以上、arithmetic の場合(Comm(Γ) は
  PSL2(R)⁰ で稠密——離散ではない、原文 §2 の対偶側)にも `Fuchsian` 型の
  値を返さねばならず、これは**論理的に無理**——`Comm` は「非arithmeticの
  場合にのみ意味を持つ部分関数」として書き直す必要がある(あるいは
  `MargulisArithmetic`/`ShimuraArithmetic` を経由した場合分けを型に
  持ち込む)。★★これは**Interface の実装可能性そのものに関わる欠陥**
  であり、単なる「橋渡しフィールドが足りない」より深い——`Comm` を
  正しく全域化するには **Margulis の二分法自体を証明する(あるいは
  posit した axiom として明示的に受け取る)** ことが避けられない。
- **`Proposition 3.2` 自身の結論**(`Γ_Z ⊆ Comm(Γ_X)`)は実は
  `Comm(Γ_X)` が離散であることを要求しない(単なる部分群包含の主張)
  ——`corrhyp_prop_3_2` が離散性抜きで証明できたのはこのため。
  ★つまり「`Proposition 3.2` を Skeleton の型のまま Found にする」ことは、
  `Comm` の全域化問題という**§2 由来の別の欠落**を先に解決しない限り、
  型すら整合しない。

**結論(2026-09-04、ユーザーとの対話を経て実施)**: 原文 p.5 の
`Comm(Γ) := {γ ∈ PSL₂(ℝ)⁰ | (γ·Γ·γ⁻¹) ∼ Γ}` を読み直すと、`Comm` の
コドメインを `Fuchsian`(離散)にしたこと自体が誤りで、原文は最初から
「`PSL₂(ℝ)⁰` の部分群」(離散性不問)として定義している。★★**Interface を
実際に修正した**(最小限、加算のみ): `Fuchsian`(離散性を要求しない)は
そのまま残し、`IsDiscrete : Fuchsian → Prop` + `Gamma_isDiscrete : ∀ X,
IsDiscrete (Gamma X)` を新設——`Comm`/`Sub`/`FiniteIndexIn` の型も
`Section2.lean`/`Section3.lean` の記述も一切変えていない(純追加、影響は
ゼロ)。

★★★**さらに一歩進めて、`HyperbolicCurveData` の具体的な項
`corrHypInstance` を `Found/CorrHyp/Instance.lean` に構成した**——
`Space := FuchsianGroup`、`Fuchsian := Subgroup(SL(2,ℝ))`(離散性不問)、
`Gamma`/`Comm`/`Sub`/`FiniteIndexIn`/`IsDiscrete`/`Gamma_isDiscrete` は
`FuchsianGroup.lean` の本物の群論で埋め、`Proposition 3.2` が読まない
残りの field(`MargulisArithmetic` 等)は型合わせの placeholder と明記した。
この `corrHypInstance` において
**`Skeleton.CorrHyp.prop_3_2` の `sorry` を文字通り埋める**
(`prop_3_2_at_instance`、`funext + rfl` で Skeleton の文と関数として
完全一致することを確認済み)ことに成功した——「関連する具体モデルの結果」
ではなく、**Skeleton の主張そのものの実装**という、これまでで最も厳密な
Found 完成(第26件)。

★★★★**2026-09-04続報(第27件)、上の候補を実際に完成させた**。
`Found/CorrHyp/ModularExample.lean`: モジュラー群 `SL(2,ℤ)`(`SL(2,ℝ)` への
像、`discreteSpecialLinearGroupIntRangeSL` で離散)と主合同部分群 `Γ(2)`
(`CongruenceSubgroup.Gamma 2`、mathlib に `FiniteIndex` 済み)を
`corrHypInstance` の `FuchsianGroup` として具体化し、`T=[[1,1],[0,1]]`
witness で両者が**相異なる**ことを確認した上で
`isIsogenous_witness : IsIsogenous corrHypInstance FG_SL2Z FG_Gamma2`
(`sorry` 無し、標準3公理のみ)を得た——`Definition 1.5` が
`corrHypInstance` において**教科書的に意味のある**非空虚性を持つことの
具体的な証拠(有限指数の伝播は `Subgroup.relIndex_map_map` で
`ker(単射) = ⊥` を消し込むだけで閉じた)。

★§2 `Definition 2.1`(`InfinitelyManyCorr`)を `Γ_SL2Z`(モジュラー群)で
非空虚に示せないか検討した(2026-09-04): `SL(2,ℤ)` は古典的に arithmetic
で、その commensurator は `PGL(2,ℚ)`(`Γ_SL2Z` の中で無限指数)——これが
示せれば `Definition 2.1` も非空虚性witnessを得られる。だが必要な
「`diag(2,1/2) ∈ Comm(Γ_SL2Z)`」(Hecke型の共役、`gΓg⁻¹ ∩ Γ` が両側
有限指数)は mathlib に対応する補題が無く(`Hecke`/対角共役で `grep` 0件)、
`Γ(2)` の例のように既存の `FiniteIndex` インスタンスを流用できない——
一から行列計算で `gΓg⁻¹ ∩ Γ` を合同部分群型の集合と同一視する必要がある。

★★★★★**2026-09-04続報(第30件)、完成した**。`g:=diag(2,1/2)`
に対し `g⁻¹Mg = [[a,b/4],[4c,d]]`・`gMg⁻¹ = [[a,4b],[c/4,d]]` を成分計算で
確認し(`conj_inv_formula`/`conj_formula`)、`Γ(4)`(mathlib に
`FiniteIndex` 済み)が両方向の有限指数下界になることを示して
`g2 ∈ Comm(Γ_SL2Z)` を確立、さらに `g2^k`(`k≥1`)の `(1,1)` 成分
`2⁻ᵏ ∈ (0,1)` が整数になりえないことから `g2^k ∉ Γ_SL2Z` を示し、
`{g2^n}` が可算無限個の相異なる剰余類を与えることで
**`Γ_SL2Z` が `Comm(Γ_SL2Z)` の中で無限指数**であることを確立した
(`Found/CorrHyp/HeckeExample.lean`、`infinitelyManyCorr_witness`、
sorry無し・標準3公理のみ)——`Definition 2.1` の非空虚性、§2 の初項目。

★配管の詰め(前回詰まっていた `Matrix.SpecialLinearGroup` の
instances-transparency 問題)は「`set A := (M : Matrix ...)` で台の
`Matrix` を先に1回だけ取り出し、以後 `A` だけで計算する」という
具体的な回避策で解決した——`tools/lean-idioms.md` に追記済み。

★★★★★★**2026-09-04続報(第31-32件)、`Definition 2.3`
(Shimura-arithmetic)を `Γ_SL2Z` について完全に実現した**。
`Found/CorrHyp/ShimuraArithmeticData.lean`: `F:=ℚ`(無限素点1個、
「他で非自明」が空虚に真)・`A:=M_2(ℚ)`(`Algebra.IsCentral.matrix`・
`IsSimpleRing.matrix`・`Module.finrank_matrix` で「4次元中心的単純多元環」
=四元数環という特徴づけを直接満たす、mathlib の基底つき
`QuaternionAlgebra` を経由しない逸脱)・`matrixEquivTensor` で
「唯一の無限素点で自明」・`O_A:=M_2(ℤ)`(`Γ_SL2Z = O_A ∩ SL_2(ℝ)` を
成分レベルで証明)という4データすべてを構成し、`Γ_SL2Z` が
Shimura-arithmetic であることを丸ごと証明した(sorry無し、標準3公理)。
`Definition 2.2`(Margulis-arithmetic)は代数群の部分群スキーム分類・
実点のリー群構造という mathlib に丸ごと不在の理論を要し、引き続き
人年規模。

★次の一手(未着手): `corrHypInstance.ShimuraArithmetic`(現状
`fun _ ↦ False` の placeholder)を `ShimuraArithmeticWitness` へ差し替え
れば、`Definition 2.3` の `.src` を正当に主張できる(道具の集計を
6/24→7/24 に動かせる見込み)。`Instance.lean`(先)→`ModularExample.lean`
→`ShimuraArithmeticData.lean`(後)というimport順のため、差し替えには
ファイル構成の見直しが要る——`prop_3_2_at_instance` 等が
`MargulisArithmetic`/`ShimuraArithmetic` の値を使わないことは確認済み
なので安全なはずだが、未検証。

`Definition 3.1`(`hyperbolicCore`)は今の `core := id` placeholder では
実装にならない(本物の商構成、`Lemma 5.1` の内容)。`Theorem 3.3` は
引き続き Riemann–Hurwitz/Gauss–Bonnet(mathlib 不在)待ち。

★★★**重要な歯止め(2026-09-04 確認)**: `corrHypInstance` の
`MargulisArithmetic`/`ShimuraArithmetic := fun _ ↦ False` という
placeholder を使って `Proposition 2.4`(`MargulisArithmetic Γ ↔
ShimuraArithmetic Γ`)や `Theorem 2.5`(`Arithmetic ↔ InfinitelyManyCorr`)
を「閉じる」ことは**してはいけない**——`Proposition 2.4` は `False ↔ False`
に退化して形式的には閉じるが、これは Margulis/Shimura-arithmeticity の
**実際の内容について何も示していない**(`corrhyp_prop_3_2` のケースと違い、
この2項目は arithmeticity の値そのものが結論に効くので、`prop_3_2` のときの
「未使用の仮定だから placeholder で害はない」という論法が通用しない——
実際 `Theorem 2.5` は my instance では**偽になる**: `SL(2,ℤ)` は現実には
arithmetic で `Comm` の中で無限指数のはずだが、`corrHypInstance` では
`Arithmetic := False` と決め打っているため矛盾する)。★これは
[[report-progress-not-shortfall]] や PGC セッションが記録した「自由な
データによる退化」と同型の罠——`.src` を付けてよいのは、placeholder が
**その項目の結論に一切影響しないことを確認した場合のみ**。

関連: [[leaf-first-with-graph-feedback]] / [[leaves-are-measured-not-guessed]] /
[[measure-mathlib-before-skeleton]] / [[genell-track-b]] / [[corrhyp-track-goal]] /
[[stale-status-read-lean-first]] / [[no-wall-decompose-instead]]

## 4. §4(`Lemma 4.1`・`Theorem 4.2`)の設計案

★§4 は §2/§3/§5/§6 と性質が違う——**mathlib に無い外部定理は不要**
(`isLimit_specKCone`、`Spec K = lim Spec R` が既に sorry 無しで完成しており、
`AffineTransitionLimit.lean` の spreading-out 定理を直接呼べる)。

### 2026-09-04 前半: 関数体案を検討して頓挫(記録として残す)

「`Space := (k 上の関数体) ⊕ (K 上の関数体)`」という和型を検討したが、
`Ext`(係数拡大)がテンソル積 `F ⊗_k K` に対応し、これが**一般には体に
ならない**ことで頓挫した(`Space` 型を全称的に保つのが重くなる)。

### 2026-09-04 後半: Scheme 案に転換、`FEt`/`Ext` を実装完了 ★★★

ユーザーから「工数が多いことで踏みとどまる必要はありません」との明示的な
指示を受け、当初「有限エタール射・ファイバー積の圏論的な整備一式が要る」と
見積もっていた **Scheme 案**に実際に着手したところ、mathlib の基盤が
見積もりよりずっと厚く、**1セッションで完成した**:

- `Found/CorrHyp/SchemeFEt.lean`: `Space := Over BaseK`(`BaseK := Spec ℚ`)。
  `FEt`(有限エタール射)は `AlgebraicGeometry.Etale`/`IsFinite`という
  mathlib の `MorphismProperty` 基盤(合成・恒等・base change 安定性が
  `instance` として自動)にそのまま乗り、`idFEt`/`comp`/`pullback`/`pbFst`/
  `pbSnd` が `infer_instance` と `MorphismProperty.pullback_fst`/`pullback_snd`
  だけで完成。
- `Ext`(係数拡大 `(-)_K`、`K := ℝ`)は、手作りの `Limits.pullback.map` ではなく
  **mathlib 自身の `Over.pullback f ⋙ Over.map f`**(`f : Spec ℝ ⟶ Spec ℚ`)
  という関手の合成として構成すると、`CategoryTheory.MorphismProperty.
  overPullbackMap` が base change 安定性からの遺伝を一発で与えてくれた
  ——`extFEt`(`Ext` の射側の作用)がこれで完成。
- `HyperbolicCurveData.Space : Type` を `Type u`(universe 多相)に変更
  (逸脱、`HyperbolicCurve.lean` に記録)——`Over BaseK : Type 1` が
  旧来の `Type`(`= Type 0`)固定と衝突したため。既存の `corrHypInstance`/
  `corrHypInstance2`(`u=0`)には影響なし。
- `Found/CorrHyp/Instance3.lean`: 上記を全部使い `corrHypInstance3 :
  HyperbolicCurveData` を具体的に構成(`Lemma 4.1` が読まないフィールドは
  安全な placeholder)。`Lemma 4.1` の statement が `corrHypInstance3` に
  対して型検査を通ることを確認。

**結論: `Space`/`FEt`/`Ext` の一貫した具体化という設計上のブロッカーは
解消した。** 残るのは `Lemma 4.1` **本体**(spreading-out の実際の証明)。

### `Lemma 4.1` 本体に向けた実測(2026-09-04、続き)

`toSchemeDiagram ℚ ℝ : (FgSubalgebra ℚ ℝ)ᵒᵖ ⥤ Scheme`(`FieldLimit.lean`)は
`AffineTransitionLimit.lean` が要求する前提を**両方とも自動で満たす**ことを
確認した:
```
example : IsCofiltered (FgSubalgebra ℚ ℝ)ᵒᵖ := inferInstance   -- OK
example {i j} (f : i ⟶ j) : IsAffineHom ((toSchemeDiagram ℚ ℝ).map f) := by
  infer_instance   -- OK
```
`AffineTransitionLimit.lean` には「極限の中で2つの射が一致する⟹有限段階で
一致する」型の spreading(`Scheme.exists_hom_hom_comp_eq_comp_of_
locallyOfFiniteType` 等)・「アフィン開被覆は有限段階に降りる」
(`Scheme.exists_isOpenCover_and_isAffine`)・「アフィン性は有限段階に降りる」
(`Scheme.exists_isAffine_of_isLimit`)は**すでにある**。

★★**まだ無いもの**(`lemma_4_1.needs` の1番目、EGA IV 8.8-8.10 相当の核心):
「極限上の**新しいスキーム**(有限型・有限エタールという条件つき)が、
ある有限段階の**新しいスキーム**の base change として存在する」という
**構成的な**降下——上の3つはどれも「すでにある対象・射」の spreading で、
「新しい対象を`Z_K`から`Z`へ**作る**」ものではない。組み立て方の見通し:
`exists_isOpenCover_and_isAffine` でアフィン開被覆を有限段階に降ろし→
各アフィン片を `RingHom.EssFiniteType.exists_eq_comp_ι_app_of_isColimit`
(`Algebra/Category/Ring/FinitePresentation.lean`)型の環準同型の spreading で
有限段階の環に降ろし→貼り合わせデータ(遷移射の一致)を
`exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType` で有限段階に降ろし→
`Etale`/`IsFinite` の条件自体も有限段階で成り立つことを言う(これは
`RingHom.Etale`/`FormallyEtale` の「有限生成部分環での近似」相当、
未確認)。**これは1つの既存補題では閉じず、組み立てが要る**——
`lemma_4_1.needs` の2番目(étale 剛性)・3番目(標数0接空間)は
まだ手つかず。

★★★2026-09-04さらに続報: `RingHom.Etale`(`FormallyUnramified ∧ Flat ∧
FinitePresentation` に分解可能、`Algebra.Etale.iff_flat_and_
formallyUnramified`)が余極限に沿って spreading する**汎用**定理は
mathlib のインデックスに見当たらない——`FinitePresentation` 部分は
`EssFiniteType.exists_eq_comp_ι_app_of_isColimit` を使い回せるが、
`Flat`/`FormallyUnramified` 部分は個別に確認が要る(未着手)。

★★★★**設計上の発見(2026-09-04、`corrHypInstance3` で `lemma_4_1` を
直接 `intro` して発覚)**: `Skeleton/CorrHyp/Section4.lean` の `lemma_4_1` は
結論を `∃ Z, ∃ (h : ZK = D.Ext Z), ...`(**命題的な等号 `=`**)で書いている。
一方 `HyperbolicCurveData.Iso : Space → Space → Prop` は witness(実際の
同型射)を一切持たない**裸の Prop**。`FuchsianGroup`(`Subgroup`、
extensional な等号を持つ)モデルではこの `=` は無害だったが、**`Scheme`
モデルでは一般に強すぎる**——`ZK` が `Ext X` と同型な**別の構成**
(例: 環を `ULift`/`Spec` で作り直したもの)であっても、finite étale な
correspondence 自体は(同型射も finite étale なので)問題なく存在しうる。
このとき結論の「`ZK = D.Ext Z` を満たす `Z` が存在する」は、`ZK` が
`Ext(-)` の形に**命題的に**一致する保証が無い限り**偽になりうる**——
これは EGA IV の降下理論を組んでも解決しない、**独立した形式化上のギャップ**
(数学的な困難ではなく、`Iso` が witness を持たないという Interface 設計の
issue)。

**次への示唆**: この issue の解消(`Iso` に実際の同型射データを持たせる、
または `lemma_4_1` の結論を `D.Iso` ベースに書き換える)は `Lemma 4.1` の
**構成**そのものを不要にはしない(`Z` を実際に作る EGA IV の仕事は
`=` でも `Iso` でも変わらず必要)——ただし構成できた後の「一致」の
証明が `rfl`/構成的等号ではなく同型の構成で済むようになる、という点で
後工程を楽にする。`Interface/CorrHyp/HyperbolicCurve.lean` を触る前に、
他の numbered item(`FinitelyManyUpTo` 経由で `Iso` を使う `Theorem 3.3`・
`Theorem 4.2`・`Theorem 2.6`)への影響も見る必要がある——`Iso` を
Prop のままにするか、witness 付きの構造体に変えるかは要検討。

次にここへ戻るときの出発点: 上の「構成的な降下」を`FEtK`(有限エタール、
`IsFinite ∧ Etale`)の場合に特化した1つの補題として切り出し、
`exists_isOpenCover_and_isAffine` から着手する。

### 2026-09-04さらに続報: `StandardEtalePair` の降下が完成 ★★★

上の「組み立て方の見通し」の最初の環——標準エタール表示(`f`・`g` の2多項式・
`cond`)の降下——を実際に `FieldLimit.lean` に実装し、sorry 無しで完成させた:

- `exists_fg_subalgebra_polynomial`(1変数)・`_pair`(2変数、同じ部分環に
  同時に降ろす)・`_pair_monic`(monic 性も遺伝)・`_family`(有限個の多項式の
  一般化)——多項式は係数が有限個なので、ある有限生成 `k`-部分環へ必ず
  降ろせるという初等的だが必要な部品。
- `exists_fg_subalgebra_standardEtaleCond`——`cond`(`f'*p₁+f*p₂=g^n`)自体は
  `f, g, p₁, p₂` を`_family`で同時に降ろした上で、`Polynomial.map`
  (`algebraMap R K` が単射なので単射)の単射性で等式を引き戻す。
- `exists_fg_subalgebra_standardEtalePair`——上の2つを束ね、
  `Algebra.StandardEtalePair K` 全体がある有限生成 `k`-部分環上の
  `StandardEtalePair` の像であることを示す。

**残る部品**(`Algebra.StandardEtalePresentation` まで完成させるには):
(a) `x : S`(`lift` が bijective になる元)自体の降下——`S` 自体が
「有限生成部分環上のある `S₀` の base change」であることも言う必要があり、
これは `S`(有限エタール多元環)を `R[X]/(f)` のこの `x` による同型で
"再定義" することで回避できる可能性がある(`x` を明示的に構成する必要が
無くなる)、(b) 一般の有限エタール多元環は「至る所で」標準エタールとは
限らず(`Algebra.basicOpen_subset_etaleLocus_iff_etale`)、**局所的に**
しか標準エタール表示を持たない——`Z_K` 全体のアフィン開被覆
(`exists_isOpenCover_and_isAffine`)の**各アフィン片**を、さらに
étale-locus のレベルで細分してから、この降下を1片ずつ適用し、最後に
`RingHom.EssFiniteType.exists_comp_map_eq_of_isColimit`(遷移射の一致の
降下)で貼り合わせる、という二重の被覆の組み立てが要る。

### 2026-09-04さらに続報: (a) の見通しを精緻化——「テンソル積の余極限保存」は不要

`exists_fg_subalgebra_standardEtalePair` で得た `f₀, g₀ : Polynomial R.1`
(`f₀.map = f`・`g₀.map = g`)から、**`S₀ := (R.1[X]/(f₀))` を `g₀` で局所化
したもの**を直接定義すればよい——これが「有限生成部分環上のある `S₀`」の
候補で、`x` を明示的に探す必要は無い(原文の `x`・`hasMap`・`lift_bijective`
は「`S` が既にこの形をしている」ことの確認用データであり、**こちらから
`S₀` を作る側では不要**)。

残る唯一の作業は「`S₀ ⊗_R K ≅ S`」——これは
`(R[X]/(f₀)) ⊗_R K ≅ K[X]/(f₀.map) = K[X]/(f)`(商環の base change、
`Ideal.map`/`Polynomial.quotientEquiv` 相当)と、局所化の base change
(`IsLocalization.baseChangeEquiv` 相当)を合成すればよい、**標準的な
可換環論の事実**——これは一般の「余極限の保存」(`Under.pushout` が
left adjoint であることを使う圏論的な議論、`Adjunction.leftAdjoint_
preservesColimits`)よりずっと軽い道具で済む可能性が高い。

★`Under.pushout`(`CategoryTheory.Under.pushoutIsLeftAdjoint`、
`PreservesColimit` が `infer_instance` で自動的に付くことは確認済み)による
「テンソル積は左随伴だから余極限を保存する」という一般論も**筋は良い**が、
`toRingCat k K` を `Under (CommRingCat.of k)` へ持ち上げる段で
`Under.mk`/`Under.homMk` の defeq 周りの配管が詰まった(`Over.mk` で
見た「`instances` 透明度」問題の `Under` 版、未解決のまま持ち越し)——
`corrhyp-goal.md` に記録だけして、上記のより軽い道(商環・局所化の base
change)を先に試す方が近道と判断した。

### 2026-09-04さらに続報: `StandardEtalePair.Ring` の base change が完成 ★★★★

上記の見通しどおり、`S₀ := R[X]/(f)` を局所化する代わりに mathlib 自身の
`StandardEtalePair.Ring`(`Polynomial (Polynomial R) ⧸ span {C f, X*C g-1}`、
明示的な二変数多項式商)を直接使う道で、base change の同型を
**完全に(sorry 無しで)構成した**:

- `Bivariate_equivMvPolynomial_map`——`Polynomial.Bivariate.equivMvPolynomial`
  (二変数多項式環 `≃ MvPolynomial (Fin 2)`)が係数の写像 `φ` と可換である
  ことを、両辺が環準同型であることから `Polynomial.ringHom_ext'` を2回
  (外側の変数・内側の変数)適用して示した。
- `standardEtalePair_ring_baseChange`——`equivMvPolynomialQuotient` が使う
  イデアル(`{C f, X*C g-1}` の span)の生成元の像が、`includeRight` →
  `MvPolynomial.algebraTensorAlgEquiv`(多変数多項式環の base change、
  mathlib 既存)で運ぶと、対応する `K` 側の生成元に**文字通り一致する**
  ことを `Ideal.map_span` + 上の可換性で計算した。
- `standardEtalePairRingBaseChange`——上の2つと
  `Algebra.TensorProduct.tensorQuotientEquiv`(商環の base change)・
  `Ideal.quotientEquivAlg`(イデアルの対応から商の同型を作る)を合成し、

  **`TensorProduct R K P₀.Ring ≃ₐ[K] (P₀.map (algebraMap R K)).Ring`**

  を得た。「余極限の保存」という一般論は不要で、mathlib の base change
  補題(`MvPolynomial.algebraTensorAlgEquiv`・`tensorQuotientEquiv`・
  `Ideal.quotientEquivAlg`)の組み合わせだけで閉じた。

**これで `Lemma 4.1` の構成的降下のうち、最も技術的に重かった「1つの
標準エタール表示の base change」が完全に閉じた。** 残るのは組み立てのみ:
(a) 一般の有限エタール多元環は「至る所で」標準エタールとは限らないため
(`Algebra.basicOpen_subset_etaleLocus_iff_etale`)、`Z_K` 全体のアフィン
開被覆(`Scheme.exists_isOpenCover_and_isAffine`)の各片をさらに
étale-locus のレベルで細分してから、この降下を1片ずつ適用すること、
(b) 各片の降下結果を `Scheme.exists_hom_hom_comp_eq_comp_of_
locallyOfFiniteType`(遷移射の一致の降下)で貼り合わせて `Z` 全体を
構成すること。数学的な核心(base change の構成)は完成したので、残るのは
**スキーム論の定型的な貼り合わせ作業**という段階に到達した。
