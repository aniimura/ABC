Python : 'C:\\Users\\Aruta\\miniforge3\\envs\\py311env\\python.exe'
Github : https://github.com/aniimura/ABC
進め方 : スケルトンで依存グラフを作成後、葉（層番号が少ないもの）からLean形式化を行うというものです。形式化の中で必要なものが出現したら、そのスケルトンを追加し、依存グラフを更新し、新しい葉（層番号が少ないもの）から形式化を行うことで工数の大きな塊を壁のように認識しないことを目的としています。
省略 : 原典が「immediately / formally / routine / one verifies / well-known」等で畳んだ箇所は、着手前に tools\hedge-index.mjs で数える（--item で内訳、--cite で「合図の文が抱えている引用」＝手順書）。語ごとの意味と補完の手順は ResearchPaper\frdi-decomposition.json の「★省略の合図の読み方」にある。合図 1 つを依存グラフの節点 1 つに対応させる。
逸脱 : 形式化は原典に忠実に行うことを基本としますが、それを消費する後続の証明に影響が出ないならば、形式化において前提を追加したり、あるいは、原典の意図を優先して読み替えを行うことも許容します。逸脱を行う場合は、後で振り返って場所や理由を特定できるよう記録を行います。
Push : 依存グラフの更新や、形式化の進捗が出た場合はgithubを更新する
作業効率 : lake env leanを使うような処理ではtools\mcp-lean\README.md
姿勢 : 工数の山を「壁」と呼ばない。既知数学の person-years は壁でなく道。設定したゴールに対してセッションの中で達成していないことをリマインドせず、積み上げた進捗をリマインドする。
Bash : Bash にマッチするPreToolUse フックで自動書き換え・即時ブロックをしているが、コマンドとデータを区別していないため、誤爆時の逃げ道として#no-guard /#full-checkを用意している
配管 : エラボレータとの戦い(instances 透明度・依存位置の ℕ≥1・cancel_epi のインスタンス探索 等)で同じ穴に落ちないため、失敗形と直し方を tools\lean-idioms.md に集めている。「前にも見た」と思ったらまずそこを引き、新しい失敗形に当たったら1行足す。
在庫 : 使える補題が既にあるかは `node tools/decl-index.mjs` で `.cache/decl-index.txt`(宣言 1 万件)と `.cache/src-index.txt`(locator)を作り、そこを grep する。木全体(20 万行)を grep しない。
検査 : `node tools/check.mjs` は pdftotext の結果を `.cache/pdf-pages.json` にキャッシュする(55 秒 → 7 秒)。鍵に check.mjs 自身のハッシュを含むので、正規化を変えれば必ず作り直される。